#' Correlation network calculation of microbial community data
#'
#' @param ps phyloseq Object, contains OTU tables, tax table and map table, represented sequences,phylogenetic tree.
#' @param method.scale sacleing method for microbiome data, rela, log,TMM,RLC···
#' @param N filter OTU tables by abundance.The defult, N=0, extract the top N number relative abundance of OTU.
#' @param r.threshold The defult, r.threshold=0.6, it represents the correlation that the absolute value
#'  of the correlation threshold is greater than 0.6. the value range of correlation threshold from 0 to 1.
#' @param p.threshold The defult, p.threshold=0.05, it represents significance threshold below 0.05.
#' @param  method method for Correlation calculation,method="pearson" is the default value.
#' The alternatives to be passed to cor are "spearman" and "kendall". recentliy "sparcc" added to method
#' @param R number of repeats for sparcc p value
#' @param ncpus number of cpu for sparcc calculate p value
#' @examples
#' data(ps)
#' result <- corMicro(ps = ps,N = 100,r.threshold=0.6,p.threshold=0.05,method = "pearson")
#' # extract cor matrix
#' cor = result[[1]]
#' @return list which contains OTU correlation matrix
#' @author Contact: Tao Wen \email{taowen@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn} yongxin liu \email{yxliu@@genetics.ac.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Tao Wen#, Penghao Xie#, Shengdie Yang, Guoqing Niu, Xiaoyu Liu, Zhexu Ding, Chao Xue, Yong-Xin Liu *, Qirong Shen, Jun Yuan*
#' ggClusterNet: an R package for microbiome network analysis and modularity-based multiple network layouts
#' iMeta 2022,DOI: \url{doi: 10.1002/imt2.32}
#' @export

corMicro = function(ps = ps,N = 0,r.threshold=0.6,
                    method.scale = "rela",
                    p.threshold=0.05,method = "pearson",R = 10,ncpus = 1){



  if (method %in% c("pearson","spearman","kendall")) {
    # ps_rela  = transform_sample_counts(ps, function(x) x / sum(x) )
    # ps_rela  = scale_micro(ps = ps,method = method.scale)
    ps_sub = filter_OTU_ps(ps = ps,Top = N)
    otu_table = as.data.frame(t(vegan_otu(ps_sub)))
    head(otu_table)
    #--- use corr.test function to calculate relation#--------
    occor = psych::corr.test(t(otu_table),use="pairwise",method=method,adjust="fdr",alpha=.05)
    occor.r = occor$r
    occor.p = occor$p

  }

  if (method %in% c("sparcc")) {
    ps_sub = filter_OTU_ps(ps = ps,Top = N)
    otu_table = as.data.frame(t(vegan_otu(ps_sub)))
    head(otu_table)
    result <- sparcc.micro(data = t(otu_table),R = R,ncpus = ncpus)
    occor.r = result[[1]]
    occor.p = result[[2]]
  }
  # occor.r[occor.p > p.threshold & abs(occor.r)<r.threshold] = 0
  occor.r[occor.p > p.threshold | abs(occor.r)<r.threshold] = 0
  return(list(occor.r,method,ps_sub,occor.p))

}





get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

reorder_cor_and_p <- function(cormat, pmat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  pmat <- pmat[hc$order, hc$order]
  list(r = cormat, p = pmat)
}



sparcc.micro <- function(
  data = data,
  R = 10,
  ncpus = 1

){
  spmatrix <- SpiecEasi::sparcc(data)

  tp0 <- proc.time()
  sp.boot <- SpiecEasi::sparccboot(
    data,
    R = R,
    ncpus = ncpus
  )
  tp1 <- proc.time()
  tp1 - tp0


  sp.p <- SpiecEasi::pval.sparccboot(sp.boot, sided = "both")



  cors <- sp.p$cors
  sp.p$pvals[is.na(sp.p$pvals)] = 1
  pvals <- sp.p$pvals
  sparCCpcors <- diag(0.5, nrow = dim(spmatrix$Cor)[1], ncol = dim(spmatrix$Cor)[1])
  sparCCpcors[upper.tri(sparCCpcors, diag=FALSE)] <- cors
  sparCCpcors <- sparCCpcors + t(sparCCpcors)

  sparCCpval <- diag(0.5, nrow = dim(spmatrix$Cor)[1], ncol = dim(spmatrix$Cor)[1])
  sparCCpval[upper.tri(sparCCpval, diag=FALSE)] <- pvals
  sparCCpval <- sparCCpval + t(sparCCpval)
  dim(sparCCpval)

  rownames(sparCCpcors) <- colnames(data)
  colnames(sparCCpcors) <- colnames(data)
  rownames(sparCCpval) <- colnames(data)
  colnames(sparCCpval) <- colnames(data)



  reordered_all_sparcc <- reorder_cor_and_p(sparCCpcors, sparCCpval)
  occor.r <- reordered_all_sparcc$r
  occor.p <- reordered_all_sparcc$p

  return(list(occor.r,occor.p))

}
