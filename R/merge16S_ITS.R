
#' merge phyloseq object of 16s and ITS
#'
#' @param ps16s phyloseq Object bacterial, contains OTU tables, tax table and map table, represented sequences,phylogenetic tree.
#' @param psITS phyloseq Object fungi, contains OTU tables, tax table and map table, represented sequences,phylogenetic tree.
#' @param N16s filter OTU tables by abundance.The defult, N=0.02, extract the top 0.02 relative abundance of OTU.
#' @param NITS filter OTU tables by abundance.The defult, N=0.02, extract the top 0.02 relative abundance of OTU.
#' @param scale transform to reletive abundance
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


merge16S_ITS <- function(ps16s = ps16s,
                         psITS = psITS,
                         N16s = 0.001,
                         NITS = 0.001,
                         scale = TRUE) {

  if (scale == TRUE) {
    ps16s  = phyloseq::transform_sample_counts(ps16s, function(x) x / sum(x) )
    psITS  = phyloseq::transform_sample_counts(psITS, function(x) x / sum(x) )


  }

  ps_16s = phyloseq::filter_taxa(ps16s, function(x) mean(x) > N16s, TRUE)#select OTUs according to  relative abundance
  ps_ITS = phyloseq::filter_taxa(psITS, function(x) mean(x) > NITS , TRUE)#select OTUs according to  relative abundance


  # ps_16s  = transform_sample_counts(ps16s, function(x) x / sum(x) )
  # ps_ITS  = transform_sample_counts(psITS, function(x) x / sum(x) )

  ### 合并细菌真菌ps对象
  otu_table_16s = as.data.frame(t(vegan_otu(ps_16s)))
  row.names(otu_table_16s) = paste("16s",row.names(otu_table_16s),sep = "_")
  otu_table_ITS = as.data.frame(t(vegan_otu(ps_ITS)))
  row.names(otu_table_ITS) = paste("ITS",row.names(otu_table_ITS ),sep = "_")

  ##由于这份文件细菌真菌用的不一样的样品吗名，所以要修改一下
  tax_table_16s = as.data.frame(vegan_tax(ps_16s))
  row.names(tax_table_16s) = paste("16s",row.names(tax_table_16s),sep = "_")
  tax_table_ITS = as.data.frame(vegan_tax(ps_ITS))
  row.names(tax_table_ITS) = paste("ITS",row.names(tax_table_ITS),sep = "_")

  #--添加一列用于注释不同的域
  tax_table_16s$filed = rep("16s",length(row.names(tax_table_16s)))
  tax_table_ITS$filed = rep("ITS",length(row.names(tax_table_ITS)))

  ##合并细菌和真菌的ps对象
  otu_table = rbind(otu_table_16s,otu_table_ITS)
  #--保证注释的表头一样
  # colnames(tax_table_16s) <- colnames(tax_table_ITS)
  tax_table = rbind(tax_table_16s,tax_table_ITS)
  dim(otu_table)

  #使用细菌的ps文件mapping 作为总得mapping文件
  mapping = as.data.frame(sample_data(ps_16s))
  head(mapping)
  # mapping$Group4 = "all_sample"
  # mapping$Group4 = as.factor(mapping$Group4)
  ##合并后的ps文件
  pallps <- phyloseq(otu_table(as.matrix(otu_table),taxa_are_rows = T),
                     sample_data(mapping),
                     tax_table(as.matrix(tax_table)))


  return(pallps )

}










