


# result= cor_Big_micro(ps = ps,N = 0,p.threshold = 0.05,r.threshold = 0.9,scale = FALSE)


cor_Big_micro = function(
  ps = ps,
  N = 0,
  p.threshold = 0.05,
  r.threshold = 0.9,
  scale = FALSE,
  method = "spearman"
  ){
  # library(igraph)
  # library(dplyr)
  # library(Hmisc)
  if (scale == TRUE) {
    ps = ps %>%
      ggClusterNet::scale_micro(method = "TMM")
  }
  x = ps %>%
    filter_OTU_ps(Top = N) %>%
    # scale_micro(method = "TMM") %>%
    vegan_otu() %>%
    t() %>%
    as.data.frame()
  occor<-WGCNA::corAndPvalue(t(x)/colSums(x),method = method)
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






