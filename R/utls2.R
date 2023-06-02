#' scale_microbiome data
#'
#' @title scaling microbiome data
#' @description will use many method for microbiome data scaling
#' @param ps phyloseq abject containg microbiome data
#' @param  method could be selected by rela, sampling, log,TMM,RLE,upperquartile
#'  et al
#' @examples
#' scale_micro(ps = ps,method = "rela")# rela, sampling, log,TMM,RLE,upperquartile

scale_micro <- function(ps,
                        method = "rela"
                        ){
  if (method == "rela") {
    ps1  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x,na.rm =TRUE) )

  }

  if (method == "sampling" ) {
    total = mean(sample_sums(ps));total
    standf = function(x,t = total)round(t*(x/sum(x)))
    ps1 = phyloseq::transform_sample_counts(ps,standf)
  }
  if (method == "log") {
    ps1  = phyloseq::transform_sample_counts(ps, function(x) log(1 + x) )

  }

  if (method %in% c("TMM","RLE","upperquartile")) {
    count = t(vegan_otu(ps))
    map <- as.data.frame(phyloseq::sample_data(ps))
    map$Group <- as.factor(map$Group)
    # create DGE list
    d = edgeR::DGEList(counts=count,group = map$Group)
    d$samples
    d = edgeR::calcNormFactors(d,method=method)
    # 生成实验设计矩阵
    design.mat = model.matrix(~ 0 + d$samples$group)
    colnames(design.mat)=levels(map$Group)
    d2 = edgeR::estimateGLMCommonDisp(d, design.mat)
    d2 = edgeR::estimateGLMTagwiseDisp(d2, design.mat)
    fit = edgeR::glmFit(d2, design.mat)

    otu = as.matrix(fit$fitted.values)
    head(otu)

    ps1 <- phyloseq::phyloseq(phyloseq::otu_table(as.matrix(otu),taxa_are_rows = TRUE),
                              phyloseq::tax_table(ps),
                              phyloseq::sample_data(ps)
    )
  }


  return(ps1)
}


#' filter_microbiome data
#'
#' @title filter microbiome data
#' @description filter microbiome data
#' @param ps phyloseq abject containg microbiome data
#' @param  Top number of top abundance tax or OTU
#' @examples
#' data(ps)
#' ps_sub = filter_OTU_ps(ps = ps,Top = 100)

filter_OTU_ps <- function(ps = ps,Top = NULL

){
  if (!is.null(Top)&Top != 0) {
    #
    ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
    # otu table
    otu_table = as.data.frame(t(vegan_otu(ps_rela)))
    otu_table$mean = rowMeans(otu_table)
    otu_table$ID = row.names(otu_table)

    otu_table<- dplyr::arrange(otu_table, desc(mean))
    subtab = head(otu_table,Top)
    otu_table2 = as.data.frame(t(vegan_otu(ps)))
    phyloseq::otu_table(ps) = phyloseq::otu_table(as.matrix(otu_table2[subtab$ID,]),taxa_are_rows = TRUE)
  } else if(Top == 0){
    ps = ps
  }

  return(ps)
}


#' filter_microbiome data upper
#'
#' @title filter microbiome data
#' @description filter microbiome data
#' @param ps phyloseq abject containg microbiome data
#' @param  filter number or percent of tax or OTU
#' @examples
#' data(ps)
#' ps_sub = filter_OTU_ps2(ps = ps,filter = 100)
# library(ggClusterNet)

filter_OTU_ps2 <- function(ps = ps,filter = NULL
){
  if (!is.null(filter)&filter != 0& filter > 1) {
    ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
    # otu table
    otu_table = as.data.frame(t(vegan_otu(ps_rela)))
    otu_table$mean = rowMeans(otu_table)
    otu_table$ID = row.names(otu_table)

    otu_table<- dplyr::arrange(otu_table, desc(mean))
    subtab = head(otu_table,filter)
    otu_table2 = as.data.frame(t(vegan_otu(ps)))
    phyloseq::otu_table(ps) = phyloseq::otu_table(as.matrix(otu_table2[subtab$ID,]),taxa_are_rows = TRUE)
  } else if(filter == 0){
    ps = ps
  } else if (filter>0 & filter < 1){
    tem = ps %>%
      phyloseq::transform_sample_counts(function(x) x / sum(x) ) %>%
      filter_taxa(function(x) sum(x ) > filter , TRUE) %>%
      vegan_otu() %>%
      t()
    phyloseq::otu_table(ps) =
      phyloseq::otu_table(as.matrix(otu_table2[row.names(tem),]),taxa_are_rows = TRUE)

  }

  return(ps)
}

selectlayout <- function(m,layout = "fruchtermanreingold"){
  if (layout == "fruchtermanreingold") {
    plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))

  }

  if (layout == "circle") {
    plotcord <-  data.frame(sna::gplot.layout.circle(m, NULL))
  }

  if (layout == "kamadakawai") {
    plotcord <-  data.frame(sna::gplot.layout.adj(m, NULL))
  }
  if (layout == "adj") {
    plotcord <-  data.frame(sna::gplot.layout.kamadakawai(m, NULL))
  }
  if (layout == "circrand") {
    plotcord <-  data.frame(sna::gplot.layout.circrand(m, NULL))
  }
  if (layout == "eigen") {
    plotcord <-  data.frame(sna::gplot.layout.eigen(m, NULL))
  }
  if (layout == "geodist") {
    plotcord <-  data.frame(sna::gplot.layout.geodist(m, NULL))
  }
  if (layout == "hall") {
    plotcord <-  data.frame(sna::gplot.layout.hall(m, NULL))
  }

  if (layout == "mds") {
    plotcord <-  data.frame(sna::gplot.layout.mds(m, NULL))
  }

  if (layout == "princoord") {
    plotcord <-  data.frame(sna::gplot.layout.princoord(m, NULL))
  }
  if (layout == "random") {
    plotcord <-  data.frame(sna::gplot.layout.random(m, NULL))
  }
  if (layout == "rmds") {
    plotcord <-  data.frame(sna::gplot.layout.rmds(m, NULL))
  }
  if (layout == "segeo") {
    plotcord <-  data.frame(sna::gplot.layout.segeo(m, NULL))
  }
  if (layout == "seham") {
    plotcord <-  data.frame(sna::gplot.layout.seham(m, NULL))
  }
  if (layout == "spring") {
    plotcord <-  data.frame(sna::gplot.layout.spring(m, NULL))
  }
  if (layout == "springrepulse") {
    plotcord <-  data.frame(sna::gplot.layout.springrepulse(m, NULL))
  }
  if (layout == "target") {
    plotcord <-  data.frame(sna::gplot.layout.target(m, NULL))
  }

  return(plotcord)
}




mergePs_Top_micro <- function(psm = psP,
                              Top = 10,
                              j = NULL){

  res <- psm %>%
    ggClusterNet::tax_glom_wt(ranks = j)

  otu = otu_table(res)
  tax = tax_table(res)

  if (is.null(j)) {
    j = ncol(tax)
  }

  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {

      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(res)= tax
  res2 <- res %>%
    ggClusterNet::tax_glom_wt(ranks = j)


  return(res2)
}

make_igraph = function(cor){
  corr <- cor
  # Construct Edge File
  edges <- data.frame(from = rep(row.names(corr), ncol(corr)),
                      to = rep(colnames(corr), each = nrow(corr)),
                      r = as.vector(corr)
  )
  # Extract half of the matrix, and remove the diagonal correlation (self related to yourself)
  edges <- dplyr::filter(edges, as.vector(lower.tri(corr)))
  colnames(edges)[3] = "weight"

  edges = edges %>% dplyr::filter(weight != 0)
  #---Set the sign of the edge
  # E.color <- edges$weight
  edges$direction <- ifelse(edges$weight>0, "pp",ifelse(edges$weight<0, "np","ns"))
  node = data.frame(name = unique(c(as.character(edges$from),as.character( edges$to))))
  row.names(node) = node$name
  # Output igraph object
  igraph  = igraph::graph_from_data_frame(edges, directed = FALSE, vertices = node)
  return(igraph)
}

# library(ggClusterNet)
# library(phyloseq)
# library(tidyverse)
# data(ps)
#
# tax = ps %>% vegan_tax() %>%
#   as.data.frame()
# head(tax)
#
# tax = remove_rankID(tax)
# tax_table(ps) = tax


remove_rankID = function(taxtab){
  taxtab$Kingdom = gsub("d__","",taxtab$Kingdom)
  taxtab$Kingdom = gsub("k__","",taxtab$Kingdom)
  taxtab$Phylum = gsub("p__","",taxtab$Phylum)
  taxtab$Class = gsub("c__","",taxtab$Class)
  taxtab$Order = gsub("o__","",taxtab$Order)
  taxtab$Family = gsub("f__","",taxtab$Family)
  taxtab$Genus = gsub("g__","",taxtab$Genus)
  taxtab$Species = gsub("s__","",taxtab$Species)
  return(taxtab)
}

# ps0 = change.OTU.name.ps(ps0)
change.OTU.name.ps = function(ps0 = ps){
  otu = ps0 %>% vegan_otu() %>% t() %>%
    as.data.frame()
  head(otu)
  rep = row.names(otu)
  row.names(otu) = paste("ASV_",1:nrow(otu),sep = "")

  tax = ps0 %>% vegan_tax() %>%
    as.data.frame()
  row.names(tax) = paste("ASV_",1:nrow(tax),sep = "")
  head(tax)
  library(Biostrings)

  tree = phy_tree(ps0)
  match(tree$tip.label,rep)

  tree$tip.label =  paste("ASV_",1:length(tree$tip.label),sep = "")

  ps = phyloseq(
    otu_table(as.matrix(otu),taxa_are_rows =TRUE),
    tax_table(as.matrix(tax)),
    sample_data(ps0),
    phy_tree(tree)
  )
  return(ps)
}




# library(ggClusterNet)
# data(ps)
# subset_samples.wt( ps = ps, Group = "Group",id = c("KO"))

subset_samples.wt = function(
    ps,
    Group,
    id,
    opst = F
){

  map = phyloseq::sample_data(ps)
  tem = map[,Group][[1]]
  tem2 = match(map[,Group][[1]],id) == 1
  if (opst == FALSE) {
    map1 = map[!is.na(tem2),]
  } else if(opst ==TRUE) {
    map1 = map[is.na(tem2),]
  }

  sample_data(ps) = map1
  return(ps)
}

# subset_taxa.wt(ps = ps,rank = "OTU",id = "ASV_1")

subset_taxa.wt = function(
    ps,
    rank,
    id,
    opst = F
){



  if (opst == FALSE) {
    if (rank %in% colnames(phyloseq::tax_table(ps))) {

      tax = ps %>%
        ggClusterNet::vegan_tax() %>%
        as.data.frame()

      tem = tax[,rank]
      tem2 = match(tax[,rank],id) == 1
      tax1 = tax[!is.na(tem2),]
      phyloseq::tax_table(ps) = phyloseq::tax_table(as.matrix(tax1))
    } else if(rank %in% c("ASV","OTU","Zotu","ZOTU")){
      otu = ps %>% vegan_otu() %>% t() %>% as.data.frame()
      otu = otu[id,]
      otu1 = otu[row.names(otu) != "NA",]
      phyloseq::otu_table(ps) =  phyloseq::otu_table(as.matrix(otu1),taxa_are_rows = TRUE)
    } else {
      print("No action")
    }
  } else if(opst ==TRUE) {
    # opst
    if (rank %in% colnames(phyloseq::tax_table(ps))) {

      tax = ps %>%
        ggClusterNet::vegan_tax() %>%
        as.data.frame()

      tem = tax[,rank]
      tem2 = match(tax[,rank],id) == 1
      tax1 = tax[is.na(tem2),]
      phyloseq::tax_table(ps) = phyloseq::tax_table(as.matrix(tax1))
    } else if(rank %in% c("ASV","OTU","Zotu","ZOTU")){
      otu = ps %>% vegan_otu() %>% t() %>% as.data.frame()
      # otu = otu[id,]
      # otu1 = otu[row.names(otu) != "NA",]
      otu1 = otu[!row.names(otu) %in% id,]
      phyloseq::otu_table(ps) =  phyloseq::otu_table(as.matrix(otu1),taxa_are_rows = TRUE)
    } else {
      print("No action")
    }

  }

  return(ps)
}



subset_taxa.wt2 = function(
    ps,
    rank,
    id,
    opst = F
){



  if (opst == FALSE) {

    if (is.null(ps@tax_table)) {
      otu = ps %>%
        ggClusterNet::vegan_otu() %>% t() %>%
        as.data.frame()
      tax = data.frame(row.names = row.names(otu),ID = row.names(otu) )
      tax_table(ps) = as.matrix(tax)
    }


    if (rank %in% colnames(phyloseq::tax_table(ps))) {

      tax = ps %>%
        ggClusterNet::vegan_tax() %>%
        as.data.frame()

      tem = tax[,rank]
      tem2 = match(tax[,rank],id) == 1
      tax1 = tax[!is.na(tem2),]
      phyloseq::tax_table(ps) = phyloseq::tax_table(as.matrix(tax1))
    } else if(rank %in% c("ASV","OTU","Zotu","ZOTU")){
      otu = ps %>% vegan_otu() %>% t() %>% as.data.frame()
      otu = otu[id,]
      otu1 = otu[row.names(otu) != "NA",]
      phyloseq::otu_table(ps) =  phyloseq::otu_table(as.matrix(otu1),taxa_are_rows = TRUE)
    } else {
      print("No action")
    }
  } else if(opst ==TRUE) {
    # opst
    if (rank %in% colnames(phyloseq::tax_table(ps))) {

      tax = ps %>%
        ggClusterNet::vegan_tax() %>%
        as.data.frame()

      tem = tax[,rank]
      tem2 = match(tax[,rank],id) == 1
      tax1 = tax[is.na(tem2),]
      phyloseq::tax_table(ps) = phyloseq::tax_table(as.matrix(tax1))
    } else if(rank %in% c("ASV","OTU","Zotu","ZOTU")){
      otu = ps %>% vegan_otu() %>% t() %>% as.data.frame()
      # otu = otu[id,]
      # otu1 = otu[row.names(otu) != "NA",]
      otu1 = otu[!row.names(otu) %in% id,]
      phyloseq::otu_table(ps) =  phyloseq::otu_table(as.matrix(otu1),taxa_are_rows = TRUE)
    } else {
      print("No action")
    }

  }

  return(ps)
}



change.rank.name = function(ps){
  tax = ps %>% vegan_tax() %>%
    as.data.frame()
  colnames(tax) =  c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  tax_table(ps) = as.matrix(tax)
  return(ps)
}


# ps.f = random.add.ps(ps = ps)

random.add.ps = function(ps = ps,add= 6,zoom = 0.3,
                         addid = "_add"
                         ){
  map = sample_data(ps)
  id = map$Group %>% unique() %>% as.character()
  A = c()
  for (i in 1:length(id)){
    ps.t = ps %>% subset_samples.wt("Group",id[i]) %>%
      filter_taxa(function(x) sum(x ) > 0 , TRUE)
    tab = phyloseq::sample_sums(ps.t) %>% mean() %>% round()
    otu  = phyloseq::rarefy_even_depth(ps.t,sample.size = tab * zoom) %>%
      filter_taxa(function(x) sum(x ) > 0 , TRUE) %>% vegan_otu() %>% t() %>%
      as.data.frame()

    head(otu)

    if (dim(otu)[2] >= add&add != 1) {
    tab.d = otu[,1:add] %>% rownames_to_column("ID")
    } else if (add == 1){
      tab.d = otu %>% rownames_to_column("ID") %>%.[,1:(add+1)]

    }

    A[1:add] = id[i]

    if (i == 1) {
      dat = tab.d
      B = A
    } else{
      dat = dat %>% full_join(tab.d,by = "ID")
      B = c(B,A)
    }

    dat2 = dat %>%
      column_to_rownames("ID")
  }

  colnames(dat2) = paste(colnames(dat2),addid,sep = "")

  mapadd = data.frame(
    row.names = colnames(dat2),
    ID = colnames(dat2),
    Group = B
  )

  head(mapadd)
  head(dat2)

  ps.f = phyloseq(
    otu_table(as.matrix(dat2),taxa_are_rows = TRUE),
    sample_data(mapadd),
    tax_table(ps)
  )
  return(ps.f)
}


#--合并ps对象
# psf = merge_amp.dat(ps.a = ps,
#                     ps.m = ps.f,
#                     rank = "OTU")

merge_amp.dat = function(
    ps.a = ps,
    ps.m = ps.f,
    method = "rela",
    rank = "OTU"
){
  #--merge
  otu1 = ps.a %>% vegan_otu() %>% t() %>%
    as.data.frame() %>% rownames_to_column("id")
  head(otu1)

  otu2 = ps.m %>% vegan_otu() %>% t() %>%
    as.data.frame()%>% rownames_to_column("id")
  head(otu2)

  otu3 = otu1 %>% inner_join(otu2,by = "id") %>% column_to_rownames("id")
  head(otu3)

  map1 = ps.a %>% sample_data() %>% as.tibble() %>%
    select(ID,Group) %>% as.data.frame()
  map2 = ps.m%>% sample_data() %>% as.tibble() %>%
    select(ID,Group) %>% as.data.frame()
  map3 = rbind(map1,map2)
  row.names(map3) = map3$ID
  tax3 = ps.a %>% vegan_tax() %>% as.data.frame()

  ps3 = phyloseq::phyloseq(
    otu_table(as.matrix(as.matrix(otu3)),taxa_are_rows = TRUE) ,
    tax_table(as.matrix(tax3)),
    sample_data(map3)
  )
  return(ps3)
}

#--基于phylsoeq对象对读数数据取整
remove_decimal = function(ps = ps){
  re = ps %>% vegan_otu()
  # t()
  #-保留小数点位数达到取整目的
  for (i in 1:dim(re)[2]) {
    re[,i] = round(re[,i],0)
  }
  #--将空缺值填充为0
  re[is.na(re)] <- 0
  otu_table(ps) = otu_table(as.matrix(t(re)),taxa_are_rows = TRUE)
  return(ps)
}


# ps = readRDS("./ps_GC.rds")
# ps1 = rm.low.area(ps = ps,threshold = 10000)
# 用于将低于某一个阈值的丰度修改为0，代谢组开发
rm.low.area = function(
    ps = ps,
    threshold = 10000
){
  dat = ps %>% vegan_otu() %>% t()
  dim(dat)
  num = dat[dat < threshold] %>% length()
  print(num < dim(dat)[1]*dim(dat)[2]/5)
  dat[dat < threshold] = 0
  # row.names(dat)
  otu_table(ps) = otu_table(as.matrix(dat),taxa_are_rows = TRUE)
  return(ps)
}

#--ps对象更换row.names
cg.Rownm.ps = function(ps,id = NULL){

  if (is.null(id)) {
    psfin = ps
  } else{
    otu = ps %>% vegan_otu() %>% t() %>%
      as.data.frame()

    head(otu)

    tax = ps %>% vegan_tax() %>%
      as.data.frame()

    head(tax)
    tax$old.rname = row.names(tax)
    row.names(tax) = tax[[id]]
    row.names(otu) = tax[[id]]

    psfin = phyloseq(
      otu_table(as.matrix(otu),taxa_are_rows = TRUE),
      tax_table(as.matrix(tax)),
      sample_data(ps)
    )
  }

  return(psfin)
}


