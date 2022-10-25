
#' Correlation network calculation of big microbial community data
#'
#' @param ps phyloseq Object, contains OTU tables, tax table and map table, represented sequences,phylogenetic tree.
#' @param N filter OTU tables by abundance.The defult, N=0, extract the top N number relative abundance of OTU.
#' @param r.threshold The defult, r.threshold=0.6, it represents the correlation that the absolute value
#'  of the correlation threshold is greater than 0.6. the value range of correlation threshold from 0 to 1.
#' @param p.threshold The defult, p.threshold=0.05, it represents significance threshold below 0.05.
#' @param method method for Correlation calculation,method="pearson" is the default value.
#' The alternatives to be passed to cor are "spearman" and "kendall".
#' @examples
#' data(ps)
#' result <- cor_Big_micro(ps = ps,N = 0,p.threshold = 0.05,r.threshold = 0.9,scale = FALSE)
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


cor_Big_micro2 = function(
    ps = ps,
    N = 0,
    p.threshold = 0.05,
    r.threshold = 0.9,
    scale = FALSE,
    method = "spearman"
){

  if (scale == TRUE) {
    ps = ps %>%
      ggClusterNet::scale_micro(method = "TMM")
  }
  x = ps %>%
    filter_OTU_ps(Top = N) %>%
    vegan_otu() %>%
    t() %>%
    as.data.frame()
  occor<-WGCNA::corAndPvalue(t(x),method = method)
  mtadj<-multtest::mt.rawp2adjp(unlist(occor$p),proc='BH')
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(t(x)/colSums(x))[2])
  ## R value
  occor.r<-occor$cor
  diag(occor.r) <- 0
  # occor.r[occor.p > 0.05|abs(occor.r)<0.4] = 0

  occor.r[occor.p > p.threshold | abs(occor.r)< r.threshold] = 0
  occor.r[is.na(occor.r)]=0

  return(list(cor = occor.r))
}






