#' Correlation network calculation of microbial community data
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

corMicro = function(ps = ps,N = 0.02,r.threshold=0.6,p.threshold=0.05,method = "pearson"){

  #--function to extract otu table and tax table#-----
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }


  #-----Relative abundance conversion-------
  ps_rela  = transform_sample_counts(ps, function(x) x / sum(x) )
  #--根据设定的阈值筛选应该包含otu#---------
  ps_sub = filter_taxa(ps_rela, function(x) mean(x ) > N , TRUE)

  #output map table
  # design = mapping= as.data.frame(sample_data(ps_sub))
  # head(design)
  #otuput otu table
  otu_table = as.data.frame(t(vegan_otu(ps_sub)))
  head(otu_table)

  # # output tax table #-----
  # tax_table = as.data.frame(vegan_tax(ps_sub ))
  # head(tax_table)



  #--- use corr.test function to calculate relation#--------
  occor = psych::corr.test(t(otu_table),use="pairwise",method=method,adjust="fdr",alpha=.05)
  occor.r = occor$r
  occor.p = occor$p

  occor.r[occor.p > p.threshold|abs(occor.r)<r.threshold] = 0


  return(list(occor.r,method,ps_sub,occor.p))

}


