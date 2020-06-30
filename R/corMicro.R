#' Correlation network calculation of microbial community data
#'
#' @param ps phyloseq Object, could contains OTU tables, tax table and map table, represented sequences.
#' @param N number object, filter OTU tables by abundance.
#' @param r.threshold Correlation threshold
#' @param p.threshold Significance threshold
#' @param  method method for Correlation calculation,method="pearson" is the default value. The alternatives to be passed to cor are "spearman" and "kendall".
#' @examples
#' data(ps)
#' result <- corMicro(ps = ps,N = 0.02,r.threshold=0.6,p.threshold=0.05,method = "pearson")
#' # extract cor matrix
#' cor = result[[1]]
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
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
  ps_sub = filter_taxa(ps_rela, function(x) sum(x ) > N , TRUE)

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


