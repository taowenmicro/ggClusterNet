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
# ps_sub = filter_OTU_ps2(ps = ps,filter = 100)
# # library(ggClusterNet

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


  if (is.null(ps@tax_table)) {
    otu = ps %>% vegan_otu()
    tax = data.frame(row.names = colnames(otu),ID = colnames(otu),
                     group = colnames(otu))
    tax_table(ps) = tax_table(as.matrix(tax))

    }




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

# ps01 = changeSamplenames(ps = ps0,newid = "newid" )
# sample_data(ps01)
changeSamplenames = function(
    ps = ps0,
    newid = "newid"  #in sample data ,was one of the colname
){
  if (!is.null(ps@otu_table)) {
    otu = ps %>% vegan_otu() %>% t() %>%
      as.data.frame()
  }

  if (!is.null(ps@sam_data)) {
    map = ps %>% sample_data() %>% as.tibble() %>%column_to_rownames(newid ) %>%
      as.data.frame()
  }

  colnames(otu) = row.names(map)
  if (!is.null(ps@tax_table)) {
    psout = phyloseq(sample_data(map),
                     otu_table(as.matrix(otu), taxa_are_rows=TRUE),
                     tax_table(ps)
    )
  } else{
    psout = phyloseq(sample_data(map),
                     otu_table(as.matrix(otu), taxa_are_rows=TRUE))
  }
  return(psout)
}

# out.ps.data(ps = ps01,path = outpath,mark = "16S")
# out.ps.data(ps = psFun2,path = outpath,mark = "ITS")
# out.ps.data(ps = psF2,path = outpath,mark = "function.k")
# out.ps.data(ps = psG,path = outpath,mark = "compounds")
out.ps.data = function(ps = ps,
                       path = "./",
                       mark = "16S"
){
  if (!is.null(ps@otu_table)) {
    otu = ps %>% vegan_otu() %>% t() %>% as.data.frame() %>% rownames_to_column("ID.orig")
    write_delim(otu,paste0(path,mark,"otu.txt"))}

  if (!is.null(ps@tax_table)) {
    tax = ps %>% vegan_tax() %>% as.data.frame() %>% rownames_to_column("ID.orig")
    write_delim(tax,paste0(path,mark,"tax.txt"))
  }

  if (!is.null(ps@sam_data)) {
    map = ps %>% sample_data() %>% as.data.frame() %>%
      rownames_to_column("ID.seq") %>% as.tibble()
    write_delim(map,paste0(path,mark,"map.txt"))
  }
}


# tax_glom_wt.upper(ps = psG,ranks = "Metabolite")
tax_glom_wt.upper <- function(ps = psG,ranks = "Metabolite") {


  if (  is.numeric(ranks)) {
    ranks <- phyloseq::rank.names(ps)[ranks]
  }
  tax = ps %>% vegan_tax() %>% as.data.frame() %>%
    mutate(conb =!!sym(ranks) )
  head(tax)
  tax_table(ps) = tax_table(as.matrix(tax))

  ranks= "conb"

  otu <- as.data.frame(t(vegan_otu(ps)))
  tax <- as.data.frame(vegan_tax(ps))

  # building group
  tax[[ranks]][is.na(tax[[ranks]])] = "Unknown"
  tax[[ranks]][tax[[ranks]] == ""] = "Unknown"
  tax[[ranks]][tax[[ranks]] == "NA"] = "Unknown"
  split <- split(otu,tax[[ranks]])
  #calculate sum by group
  apply <- lapply(split,function(x)colSums(x[]))
  # result chack
  otucon <- do.call(rbind,apply)

  taxcon <- tax[1:match(ranks,colnames(tax))]
  taxcon <- taxcon[!duplicated(tax[[ranks]]),]

  if (is.vector(taxcon)) {
    taxcon = data.frame(row.names = taxcon,ranks = taxcon)
    colnames(taxcon) = ranks
  }

  #-tax name with NA wound be repeated with unknown
  taxcon[[ranks]][is.na(taxcon[[ranks]])] = "unknown"
  row.names(taxcon) <- taxcon[[ranks]]


  row.names(taxcon) <- gsub("-","_",row.names(taxcon))
  row.names(taxcon) <- gsub("[/]","_",row.names(taxcon))
  row.names(taxcon) <- gsub("[(]","_",row.names(taxcon))
  row.names(taxcon) <- gsub("[)]","_",row.names(taxcon))
  row.names(taxcon) <- gsub("[:]","_",row.names(taxcon))
  row.names(taxcon) <- gsub("[[]","",row.names(taxcon))
  row.names(taxcon) <- gsub("[]]","_",row.names(taxcon))
  row.names(taxcon) <- gsub("[#]","_",row.names(taxcon))
  row.names(taxcon) <- gsub("[+]","_",row.names(taxcon))
  # row.names(taxcon) <- gsub("[ ]","_",row.names(taxcon))
  row.names(taxcon) <- gsub("[']","_",row.names(taxcon))
  # row.names(taxcon) <- gsub("[']","_",row.names(taxcon))


  row.names(otucon) <- gsub("-","_",row.names(otucon))
  row.names(otucon) <- gsub("[/]","_",row.names(otucon))
  row.names(otucon) <- gsub("[(]","_",row.names(otucon))
  row.names(otucon) <- gsub("[)]","_",row.names(otucon))
  row.names(otucon) <- gsub("[:]","_",row.names(otucon))
  row.names(otucon) <- gsub("[[]","",row.names(otucon))
  row.names(otucon) <- gsub("[]]","_",row.names(otucon))
  row.names(otucon) <- gsub("[#]","_",row.names(otucon))
  row.names(otucon) <- gsub("[+]","_",row.names(otucon))
  # row.names(otucon) <- gsub(" ","_",row.names(otucon))
  row.names(otucon) <- gsub("[']","_",row.names(otucon))


  pscon <- phyloseq::phyloseq(
    phyloseq::otu_table( as.matrix(otucon),taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(taxcon)),
    phyloseq::sample_data(ps)
  )
  return(pscon)
}



tax_glom_wt.3 <- function(ps = ps,ranks = "Phylum") {




  if (  is.numeric(ranks)) {
    ranks <- phyloseq::rank.names(ps)[ranks]
  }

  tax = ps %>% vegan_tax() %>% as.data.frame() %>%
    mutate(conb =!!sym(ranks) )
  head(tax)
  tax_table(ps) = tax_table(as.matrix(tax))

  ranks= "conb"

  otu <- as.data.frame(t(vegan_otu(ps)))
  tax <- as.data.frame(vegan_tax(ps))

  # building group
  tax[[ranks]][is.na(tax[[ranks]])] = "Unknown"
  tax[[ranks]][tax[[ranks]] == ""] = "Unknown"
  tax[[ranks]][tax[[ranks]] == "NA"] = "Unknown"
  split <- split(otu,tax[[ranks]])
  #calculate sum by group
  apply <- lapply(split,function(x)colSums(x[]))
  # result chack
  otucon <- do.call(rbind,apply)

  taxcon <- tax[1:match(ranks,colnames(tax))]
  taxcon <- taxcon[!duplicated(tax[[ranks]]),]

  if (is.vector(taxcon)) {
    taxcon = data.frame(row.names = taxcon,ranks = taxcon)
    colnames(taxcon) = ranks
  }

  #-tax name with NA wound be repeated with unknown
  taxcon[[ranks]][is.na(taxcon[[ranks]])] = "unknown"
  row.names(taxcon) <- taxcon[[ranks]]


  pscon <- phyloseq::phyloseq(
    phyloseq::otu_table( as.matrix(otucon),taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(taxcon)),
    phyloseq::sample_data(ps)
  )
  return(pscon)
}


# library(ggClusterNet)
# library(tidyverse)
# library(phyloseq)

scale_omics = function(
    ps = ps,
    method = "TPM",
    efflength =  1000

){

  otu = ps %>% vegan_otu()
  if (method == "TPM") {
    efflength =  efflength
  }

  if (method == "TPM") {
    counts2TPM <- function(count=count, efflength=  efflength){
      RPK <- count/(efflength/1000)   #每千碱基reads (Reads Per Kilobase) 长度标准化
      PMSC_rpk <- sum(RPK)/1e6        #RPK的每百万缩放因子 (“per million” scaling factor ) 深度标准化
      RPK/PMSC_rpk
    }

    tpm <- as.data.frame(apply(otu,2,counts2TPM) %>% t())
    otu_table(ps) = otu_table(as.matrix(tpm),taxa_are_rows = TRUE)
  }

  if (method == "FPKM") {
    expr1 = otu/efflength
    fpkm = t(t(expr1)/colSums(otu)) * 10^9
    otu_table(ps) = otu_table(as.matrix(fpkm),taxa_are_rows = TRUE)
  }
  if (method == "RPKM") {
    RPKM = function(count=count, efflength=  efflength){
      lib.size = sum(count)
      rpm = (count*10^9)/lib.size
      rpkm = rpm/efflength
      return(rpkm)
    }

    rpkm = RPKM(count= otu, efflength=  efflength) %>% t()

    otu_table(ps) = otu_table(as.matrix(rpkm),taxa_are_rows = TRUE)

  }

  return(ps)
}



rowSD = function(x){
  apply(x,1, sd)
}

rowCV = function(x){
  rowSD(x)/rowMeans(x)
}


add.id.facet = function(x = d,group = "id"){
  tem = x %>% distinct(!!sym(group))
  colnames(tem) = "nm.tem"
  id.s = tem$nm.tem
  for (i in 1:length(id.s)) {
    tem2 = x %>% filter(!!sym(group) == id.s[i])
    tem2$id.facet = paste0(id.s[i],"_",1:dim(tem2)[1])

    if (i == 1) {
      tem3 = tem2
    } else{
      tem3 = rbind(tem3,tem2)
    }
  }
  return(tem3)
}



line.across.facets.network0 <- function(p, from=1, to=2,
                                       from_point_id=1,
                                       to_point_id=1,
                                       plotout = F,
                                       gp=gpar(lty=2, alpha=0.5)){
  if (TRUE %in% grepl("ggplot", class(p))) {
    g <- ggplot_gtable(ggplot_build(p))
  } else {
    g <- p
  }

  # if one of them contains NA, return g
  if (NA %in% c(from, to, from_point_id, to_point_id)) return(g)

  # collect panel viewport names and index numbers in the grob
  panel_vps <- c()
  id_n <- c()
  for (i in 1:length(g$grobs)) {
    if (str_detect(g$layout[i, "name"], "panel") & g$grobs[[i]]$name != "NULL") {
      p_name <- g$layout[i, "name"]
      panel_vps <- c(panel_vps, p_name)
      id_n <- c(id_n, i)
    }
  }

  panel_vps %>%
    str_replace("panel-", "") %>%
    str_split("[\\-\\.]") %>%
    map_chr(1) -> ind_col
  ind_col <- as.numeric(ind_col)

  panel_vps %>%
    str_replace("panel-", "") %>%
    str_split("[\\-\\.]") %>%
    map_chr(2) -> ind_row
  ind_row <- as.numeric(ind_row)

  my_dim <- c(max(ind_row), max(ind_col))
  x <- 1:length(id_n)
  L <- length(x)
  # x[(L+1):(my_dim[1]*my_dim[2])] <- NA
  m1 <- as.vector(matrix(x, nrow=my_dim[1], byrow=TRUE))

  x2 <- 1:L
  xx <- as.vector(!is.na(m1))
  xx[xx] <- x2
  xx[!xx] <- NA
  m2 <- as.vector(matrix(xx, nrow=my_dim[1]))

  from <- m2[m1==from]
  from <- from[complete.cases(from)]
  to <- m2[m1==to]
  to <- to[complete.cases(to)]

  # select points to be connected
  pnames1 <- names(g$grobs[[id_n[from]]]$children)
  pnames2 <- names(g$grobs[[id_n[to]]]$children)

  pname1 <- pnames1[str_detect(pnames1, "geom_point.points")]
  pname2 <- pnames2[str_detect(pnames2, "geom_point.points")]

  p1 <- g$grobs[[id_n[from]]]$children[[pname1[1]]]
  p2 <- g$grobs[[id_n[to]]]$children[[pname2[1]]]

  g1 <- with(g$layout[id_n[from],],
             gtable_add_grob(g,
                             moveToGrob(p1$x[from_point_id],
                                        p1$y[from_point_id]),
                             t=t, l=l,name = paste0(from,to,from_point_id,to_point_id)))
  g2 <- with(g1$layout[id_n[to],],
             gtable_add_grob(g1,
                             lineToGrob(p2$x[to_point_id],
                                        p2$y[to_point_id], gp=gp),
                             t=t, l=l,name = paste0(from,to,from_point_id,to_point_id)))

  # print(p1$x[from_point_id])
  # print(p1$y[from_point_id])
  # print(p2$x[to_point_id])
  # print(p2$y[to_point_id])
  # g <- gtable_add_grob(g,
  #                      moveToGrob(p1$x[from_point_id],
  #                                 p1$y[from_point_id]),
  #                      t=id_n[from], l=id_n[from])
  #
  # g <- gtable_add_grob(g,
  #                      lineToGrob(p2$x[to_point_id],
  #                                 p2$y[to_point_id], gp=gp),
  #                      t=id_n[2], l=id_n[2])
  # tried curve, but it seems no function for curve across viewports, this is within viewport plot
  # g <- with(g$layout[id_n[to],],
  #           gtable_add_grob(g,
  #                           segmentsGrob(p1$x[from_point_id],
  #                                     p1$y[from_point_id],
  #                                     p2$x[to_point_id],
  #                                     p2$y[to_point_id], gp=gp),
  #                           t=t, l=l))


  g2$layout$clip <- "off"
  if (plotout==TRUE) grid.draw(g2)
  return(g2)
}
line.across.facets.network <- function(p, from=1, to=2,
                                       from_point_id=1,
                                       to_point_id=1,
                                       plotout = F,
                                       gp=gpar(lty=2, alpha=0.5)){
  if (TRUE %in% grepl("ggplot", class(p))) {
    g <- ggplot_gtable(ggplot_build(p))
  } else {
    g <- p
  }

  # if one of them contains NA, return g
  if (NA %in% c(from, to, from_point_id, to_point_id)) return(g)

  # collect panel viewport names and index numbers in the grob
  panel_vps <- c()
  id_n <- c()
  for (i in 1:length(g$grobs)) {
    if (str_detect(g$layout[i, "name"], "panel") & g$grobs[[i]]$name != "NULL") {
      p_name <- g$layout[i, "name"]
      panel_vps <- c(panel_vps, p_name)
      id_n <- c(id_n, i)
    }
  }

  panel_vps %>%
    str_replace("panel-", "") %>%
    str_split("[\\-\\.]") %>%
    map_chr(1) -> ind_col
  ind_col <- as.numeric(ind_col)

  panel_vps %>%
    str_replace("panel-", "") %>%
    str_split("[\\-\\.]") %>%
    map_chr(2) -> ind_row
  ind_row <- as.numeric(ind_row)

  my_dim <- c(max(ind_row), max(ind_col))
  x <- 1:length(id_n)
  L <- length(x)
  # x[(L+1):(my_dim[1]*my_dim[2])] <- NA
  m1 <- as.vector(matrix(x, nrow=my_dim[1], byrow=TRUE))

  x2 <- 1:L
  xx <- as.vector(!is.na(m1))
  xx[xx] <- x2
  xx[!xx] <- NA
  m2 <- as.vector(matrix(xx, nrow=my_dim[1]))

  from <- m2[m1==from] %>% unique()
  from <- from[complete.cases(from)]
  to <- m2[m1==to]
  to <- to[complete.cases(to)]

  # select points to be connected
  pnames1 <- names(g$grobs[[id_n[from]]]$children)
  pnames2 <- names(g$grobs[[id_n[to]]]$children)

  pname1 <- pnames1[str_detect(pnames1, "geom_point.points")]
  pname2 <- pnames2[str_detect(pnames2, "geom_point.points")]

  p1 <- g$grobs[[id_n[from]]]$children[[pname1[1]]]
  p2 <- g$grobs[[id_n[to]]]$children[[pname2[1]]]

  g1 <- with(g$layout[id_n[from],],
             gtable_add_grob(g,
                             moveToGrob(p1$x[from_point_id],
                                        p1$y[from_point_id]),
                             t=t, l=l,name = paste0(from,to,from_point_id,to_point_id)))
  g2 <- with(g1$layout[id_n[to],],
             gtable_add_grob(g1,
                             lineToGrob(p2$x[to_point_id],
                                        p2$y[to_point_id], gp=gp),
                             t=t, l=l,name = paste0(from,to,from_point_id,to_point_id)))

  # print(p1$x[from_point_id])
  # print(p1$y[from_point_id])
  # print(p2$x[to_point_id])
  # print(p2$y[to_point_id])
  # g <- gtable_add_grob(g,
  #                      moveToGrob(p1$x[from_point_id],
  #                                 p1$y[from_point_id]),
  #                      t=id_n[from], l=id_n[from])
  #
  # g <- gtable_add_grob(g,
  #                      lineToGrob(p2$x[to_point_id],
  #                                 p2$y[to_point_id], gp=gp),
  #                      t=id_n[2], l=id_n[2])
  # tried curve, but it seems no function for curve across viewports, this is within viewport plot
  # g <- with(g$layout[id_n[to],],
  #           gtable_add_grob(g,
  #                           segmentsGrob(p1$x[from_point_id],
  #                                     p1$y[from_point_id],
  #                                     p2$x[to_point_id],
  #                                     p2$y[to_point_id], gp=gp),
  #                           t=t, l=l))


  g2$layout$clip <- "off"
  if (plotout==TRUE) grid.draw(g2)
  return(g2)
}


remove.zero = function(ps) {
  pst = ps %>%
    # scale_micro() %>%
    filter_taxa(function(x) sum(x ) > 0 , TRUE)
  return(pst)
}

ps.group.mean = function(ps,group = "Group" ){
  otu_table = as.data.frame(t(vegan_otu(ps)))
  head(otu_table)

  design = as.data.frame(sample_data(ps))
  ## 计算相对丰度，计算每个物种丰度均值，按照均值排序
  OTU = as.matrix(otu_table)
  norm = t(t(OTU)/colSums(OTU,na=TRUE)) #* 100 # normalization to total 100
  norma = norm %>%
    t() %>% as.data.frame()
  #数据分组计算平均值
  iris.split <- split(norma,as.factor(design$Group))

  iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
  # 组合结果
  norm2 <- do.call(rbind,iris.apply)%>%
    t()
  norm2 = as.data.frame(norm2)
  norm2$mean=apply(norm2,1,mean)
  norm2$ID = row.names(norm2)
  return(norm2)
}

line.across.facets.network2 <- function(p, from=1, to=2,
                                        from_point_id=1,
                                        to_point_id=1,
                                        plotout = F,
                                        gp=gpar(lty=2, alpha=0.5)){
  if (TRUE %in% grepl("ggplot", class(p))) {
    g <- ggplot_gtable(ggplot_build(p))
  } else {
    g <- p
  }

  # if one of them contains NA, return g
  if (NA %in% c(from, to, from_point_id, to_point_id)) return(g)

  # collect panel viewport names and index numbers in the grob
  panel_vps <- c()
  id_n <- c()
  for (i in 1:length(g$grobs)) {
    if (str_detect(g$layout[i, "name"], "panel") & g$grobs[[i]]$name != "NULL") {
      p_name <- g$layout[i, "name"]
      panel_vps <- c(panel_vps, p_name)
      id_n <- c(id_n, i)
    }
  }

  panel_vps %>%
    str_replace("panel-", "") %>%
    str_split("[\\-\\.]") %>%
    map_chr(1) -> ind_col
  ind_col <- as.numeric(ind_col)

  panel_vps %>%
    str_replace("panel-", "") %>%
    str_split("[\\-\\.]") %>%
    map_chr(2) -> ind_row
  ind_row <- as.numeric(ind_row)

  my_dim <- c(max(ind_row), max(ind_col))
  x <- 1:length(id_n)
  L <- length(x)
  # x[(L+1):(my_dim[1]*my_dim[2])] <- NA
  m1 <- as.vector(matrix(x, nrow=my_dim[1], byrow=TRUE))

  x2 <- 1:L
  xx <- as.vector(!is.na(m1))
  xx[xx] <- x2
  xx[!xx] <- NA
  m2 <- as.vector(matrix(xx, nrow=my_dim[1]))

  from <- m2[m1==from] %>% unique()
  from <- from[complete.cases(from)]
  to <- m2[m1==to]
  to <- to[complete.cases(to)]

  # select points to be connected
  pnames1 <- names(g$grobs[[id_n[from]]]$children)
  pnames2 <- names(g$grobs[[id_n[to]]]$children)

  pname1 <- pnames1[str_detect(pnames1, "geom_rect.rect")]
  pname2 <- pnames2[str_detect(pnames2, "geom_rect.rect")]

  p1 <- g$grobs[[id_n[from]]]$children[[pname1[1]]]
  p2 <- g$grobs[[id_n[to]]]$children[[pname2[1]]]

  g1 <- with(g$layout[id_n[from],],
             gtable_add_grob(g,
                             moveToGrob(p1$x[from_point_id],
                                        p1$y[from_point_id]),
                             t=t, l=l,name = paste0(from,to,from_point_id,to_point_id)))
  g2 <- with(g1$layout[id_n[to],],
             gtable_add_grob(g1,
                             lineToGrob(p2$x[to_point_id],
                                        p2$y[to_point_id], gp=gp),
                             t=t, l=l,name = paste0(from,to,from_point_id,to_point_id)))

  # print(p1$x[from_point_id])
  # print(p1$y[from_point_id])
  # print(p2$x[to_point_id])
  # print(p2$y[to_point_id])
  # g <- gtable_add_grob(g,
  #                      moveToGrob(p1$x[from_point_id],
  #                                 p1$y[from_point_id]),
  #                      t=id_n[from], l=id_n[from])
  #
  # g <- gtable_add_grob(g,
  #                      lineToGrob(p2$x[to_point_id],
  #                                 p2$y[to_point_id], gp=gp),
  #                      t=id_n[2], l=id_n[2])
  # tried curve, but it seems no function for curve across viewports, this is within viewport plot
  # g <- with(g$layout[id_n[to],],
  #           gtable_add_grob(g,
  #                           segmentsGrob(p1$x[from_point_id],
  #                                     p1$y[from_point_id],
  #                                     p2$x[to_point_id],
  #                                     p2$y[to_point_id], gp=gp),
  #                           t=t, l=l))


  g2$layout$clip <- "off"
  if (plotout==TRUE) grid.draw(g2)
  return(g2)
}



merge.ps <- function(ps1 ,
                     ps2,
                     N1 = 100,
                     N2 = 100,
                     scale = TRUE,
                     onlygroup = FALSE,#不进行列合并，只用于区分不同域
                     dat1.lab = "bac",
                     dat2.lab = "fun") {

  if (scale == TRUE) {
    if (!is.null(ps16s)) {
      ps1  = phyloseq::transform_sample_counts(ps1, function(x) x / sum(x) )
    }
    if (!is.null(psITS)) {
      ps2  = phyloseq::transform_sample_counts(ps2, function(x) x / sum(x) )
    }
  }
  if (!is.null(ps1)) {
    # ps_16s = phyloseq::filter_taxa(ps16s, function(x) mean(x) > N16s, TRUE)#select OTUs according to  relative abundance
    ps_16s  =  filter_OTU_ps(ps = ps1,Top = N1)
    ###
    otu_table_16s = as.data.frame(t(vegan_otu(ps_16s)))

    if (is.null(dat1.lab)) {
      row.names(otu_table_16s) = row.names(otu_table_16s)
    } else{
      row.names(otu_table_16s) = paste(dat1.lab,row.names(otu_table_16s),sep = "_")
    }


    ## change the OTU name of bac and fungi OTU table
    tax_table_16s = as.data.frame(vegan_tax(ps_16s))
    #-- add a col marked the bac and fungi
    if ("filed" %in%colnames(tax_table_16s)) {

    } else{


      if (is.null(dat1.lab)) {
        row.names(tax_table_16s) = row.names(tax_table_16s)
      } else{
        row.names(tax_table_16s) = paste(dat1.lab,row.names(tax_table_16s),sep = "_")
      }

      tax_table_16s$filed = rep(dat1.lab,length(row.names(tax_table_16s)))

    }

  }
  if (!is.null(ps2)) {
    # ps_ITS = phyloseq::filter_taxa(psITS, function(x) mean(x) > NITS , TRUE)#select OTUs according to  relative abundance
    ps_ITS = filter_OTU_ps(ps = ps2,Top = N2)
    otu_table_ITS = as.data.frame(t(vegan_otu(ps_ITS)))
    row.names(otu_table_ITS) = paste(dat2.lab,row.names(otu_table_ITS ),sep = "_")
    tax_table_ITS = as.data.frame(vegan_tax(ps_ITS))

    if (is.null(dat2.lab)) {
      row.names(tax_table_ITS) = row.names(tax_table_ITS)
    } else{
      row.names(tax_table_ITS) = paste(dat2.lab,row.names(tax_table_ITS),sep = "_")
    }



    row.names(tax_table_ITS) = paste(dat2.lab,row.names(tax_table_ITS),sep = "_")
    tax_table_ITS$filed = rep(dat2.lab,length(row.names(tax_table_ITS)))

    if ("filed" %in%colnames(tax_table_ITS)) {
    } else{

      if (is.null(dat2.lab)) {
        row.names(tax_table_ITS) = row.names(tax_table_ITS)
      } else{
        row.names(tax_table_ITS) = paste(dat2.lab,row.names(tax_table_ITS),sep = "_")
      }

      # row.names(tax_table_ITS) = paste(dat2.lab,row.names(tax_table_ITS),sep = "_")
      tax_table_ITS$filed = rep(dat2.lab,length(row.names(tax_table_ITS)))
    }
  }


  if (!is.null(ps2) & !is.null(ps1) ) {
    ## merge OTU table of bac and fungi



    otu_table = rbind(otu_table_16s[,intersect(names(otu_table_ITS),names(otu_table_16s))],

                      otu_table_ITS[,intersect(names(otu_table_ITS),names(otu_table_16s))])

    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_16s,tax_table_ITS)
      dim(otu_table)
    } else if(onlygroup == TRUE){
      tax_table = data.frame(filed = c(tax_table_16s$filed,tax_table_ITS$filed),row.names = row.names(otu_table),id = row.names(otu_table))
    }
    #on of map table as final map table

    mapping = as.data.frame( phyloseq::sample_data(ps_16s))
    head(mapping)
    # mapping$Group4 = "all_sample"
    # mapping$Group4 = as.factor(mapping$Group4)
    ##merge all abject of phyloseq
    pallps <-  phyloseq::phyloseq( phyloseq::otu_table(as.matrix(otu_table),taxa_are_rows = TRUE),
                                   phyloseq::sample_data(mapping),
                                   phyloseq::tax_table(as.matrix(tax_table)))


  } else if(is.null(psITS) & !is.null(ps16s) ) {
    otu_table = rbind(otu_table_16s)

    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_16s)
      dim(otu_table)
    } else if(onlygroup == TRUE){
      tax_table = data.frame(filed = c(tax_table_16s$filed,row.names = row.names(otu_table),id = row.names(otu_table)))
    }
    #on of map table as final map table
    mapping = as.data.frame(sample_data(ps_16s))
    head(mapping)
    # mapping$Group4 = "all_sample"
    # mapping$Group4 = as.factor(mapping$Group4)
    ##merge all abject of phyloseq
    pallps <-  phyloseq::phyloseq( phyloseq::otu_table(as.matrix(otu_table),taxa_are_rows = TRUE),
                                   phyloseq::sample_data(mapping),
                                   phyloseq::tax_table(as.matrix(tax_table)))


  } else if (!is.null(ps2) & is.null(ps1)){
    otu_table = rbind(otu_table_ITS)

    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_ITS)
      dim(otu_table)
    } else if(onlygroup == TRUE){
      tax_table = data.frame(filed = c(tax_table_ITS$filed),row.names = row.names(otu_table),id = row.names(otu_table))
    }
    #on of map table as final map table
    mapping = as.data.frame( phyloseq::sample_data(psITS))
    head(mapping)
    # mapping$Group4 = "all_sample"
    # mapping$Group4 = as.factor(mapping$Group4)
    ##merge all abject of phyloseq
    pallps <-  phyloseq::phyloseq( phyloseq::otu_table(as.matrix(otu_table),taxa_are_rows = T),
                                   phyloseq::sample_data(mapping),
                                   phyloseq::tax_table(as.matrix(tax_table)))

  }

  tax = pallps %>% vegan_tax() %>%
    as.data.frame() %>% dplyr::select(filed,everything())
  phyloseq::tax_table(pallps) = as.matrix(tax)


  return(pallps)
}



