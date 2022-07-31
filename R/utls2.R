
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
    ps1  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )

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
# library(ggClusterNet)
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

# library(ggClusterNet)
# data(ps)
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
  otu = otu_table(psm)
  tax = tax_table(psm)

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
  phyloseq::tax_table(psm)= tax

  res <- psm %>%
    ggClusterNet::tax_glom_wt(ranks = j)

  return(res)
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
