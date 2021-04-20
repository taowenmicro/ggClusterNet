#' scale_microbiome data
#'
#' @title scaling microbiome data
#' @description will use many method for microbiome data scaling
#' @param ps phyloseq abject containg microbiome data
#' @param  method could be selected byy rela, sampling, log,TMM,RLE,upperquartile
#'  et al
#' @examples
#' scale_micro(ps = ps,method = "rela")


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
    map <- as.data.frame(sample_data(ps))
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

    ps1 <- phyloseq(otu_table(as.matrix(otu),taxa_are_rows = TRUE),
                   tax_table(ps),
                   sample_data(ps)
    )
  }


  return(ps1)
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


