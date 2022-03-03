#' Inter-Domain Ecological Network, we call this biostripe network
#'
#' @param ps phyloseq Object, contains OTU tables, tax table and map table, represented sequences,phylogenetic tree.
#' @param N filter OTU tables by abundance.The defult, N=0.02, extract the top 0.02 relative abundance of OTU.
#' @param r.threshold The defult, r.threshold=0.6, it represents the correlation that the absolute value
#'  of the correlation threshold is greater than 0.6. the value range of correlation threshold from 0 to 1.
#' @param p.threshold The defult, p.threshold=0.05, it represents significance threshold below 0.05.
#' @param  method method for Correlation calculation,method="pearson" is the default value. The alternatives to be passed to cor are "spearman" and "kendall".
#' @examples
#' data(ps)
#' result <- corMicro(ps = ps,N = 0.02,r.threshold=0.6,p.threshold=0.05,method = "pearson")
#' # extract cor matrix
#' cor = result[[1]]
#' @return list which contains OTU correlation matrix
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export




corBiostripeBig = function(data = NULL, group = NULL,ps = NULL,r.threshold=0.6,p.threshold=0.05,method = "spearman"){
  
  if (is.null(data)&is.null(group)&!is.null(ps)) {
    # otu_table = as.data.frame(t(vegan_otu(ps)))
 
    x = ps %>%
      # filter_OTU_ps(Top = N) %>%
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
    # occor.r[occor.p > p.threshold|abs(occor.r)<r.threshold] = 0
    occor.r[abs(occor.r)<r.threshold] = 0
    occor.r[occor.p > p.threshold] = 0
    tax = as.data.frame((vegan_tax(ps)))
    head(tax)
    A <- levels(as.factor(tax$filed))
    A
    # i = 1
    for (i in 1:length(A)) {
      fil <- intersect(row.names(occor.r),as.character(row.names(tax)[tax$filed == A[i]]))
      a <- row.names(occor.r) %in% fil
      occor.r[a,a] = 0
      occor.p[a,a] = 1
    }
    
  }
  
  
  if (is.null(ps)&!is.null(data)&!is.null(group)) {
    cordata <- t(data[-1])
    colnames(cordata) =data[[1]]
    #--- use corr.test function to calculate relation#--------
    # occor = psych::corr.test(cordata,use="pairwise",method=method,adjust="fdr",alpha=.05)
    # occor.r = occor$r
    # occor.p = occor$p
    
    x = cordata %>% t()
    occor<-WGCNA::corAndPvalue(t(x)/colSums(x),method = method)
    mtadj<-multtest::mt.rawp2adjp(unlist(occor$p),proc='BH')
    adpcor<-mtadj$adjp[order(mtadj$index),2]
    occor.p<-matrix(adpcor,dim(t(x)/colSums(x))[2])
    ## R value
    occor.r<-occor$cor
    
    
    #-filter--cor value
    occor.r[occor.p > p.threshold|abs(occor.r)<r.threshold] = 0
    
    #--biostripe network filter
    A <- levels(as.factor(group$Group))
    
    
    
    for (i in 1:length(A)) {
      fil <- intersect(row.names(occor.r),as.character(group[[1]][group$Group == A[i]]))
      a <- row.names(occor.r) %in% fil
      occor.r[a,a] = 0
      occor.p[a,a] = 1
    }
    
  }
  
  
  if (!is.null(ps)&!is.null(data)&!is.null(group)) {
    otu_table = as.data.frame(t(vegan_otu(ps)))
    cordata <- (data[-1])
    row.names(cordata) = data[[1]]
    
    if (!is.na(match(colnames(otu_table) , data[[1]]))) {
      cordata = t(cordata)
    }
    dim(cordata)
    dim(otu_table)
    finaldata <- rbind(otu_table,cordata)
    
    
    #--- use corr.test function to calculate relation#--------
    # occor = psych::corr.test(t(finaldata),use="pairwise",method= method ,adjust="fdr",alpha=.05)
    # occor.r = occor$r
    # occor.p = occor$p
    # occor.r[occor.p > p.threshold|abs(occor.r)<r.threshold] = 0
    
    x = finaldata
    occor<-WGCNA::corAndPvalue(t(x)/colSums(x),method = method)
    mtadj<-multtest::mt.rawp2adjp(unlist(occor$p),proc='BH')
    adpcor<-mtadj$adjp[order(mtadj$index),2]
    occor.p<-matrix(adpcor,dim(t(x)/colSums(x))[2])
    
    
    
    tax = as.data.frame((vegan_tax(ps)))
    head(tax)
    if (length(tax$filed) != 0) {
      A1 <- levels(as.factor(tax$filed))
      A1
      A2 <- levels(as.factor(group[[2]]))
      A2
      A = c(A1,A2)
      
      group2 <- data.frame(SampleID = row.names(tax),Group = tax$filed)
    } else {
      A1 = "OTU"
      A2 <- levels(as.factor(group[[2]]))
      A2
      A = c(A1,A2)
      group2 <- data.frame(SampleID = row.names(tax),Group = "OTU")
    }
    
    
    
    # i = 5
    colnames(group) = c("SampleID","Group")
    finalgru = rbind(group,group2)
    
    for (i in 1:length(A)) {
      fil <- intersect(row.names(occor.r),as.character(as.character(finalgru$SampleID)[as.character(finalgru$Group) == A[i]]))
      a <- row.names(occor.r) %in% fil
      occor.r[a,a] = 0
      occor.p[a,a] = 1
    }
    
  }
  return(list(occor.r, method, occor.p, A))
  
}











