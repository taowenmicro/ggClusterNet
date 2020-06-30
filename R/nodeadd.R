#' Use otu forms and annotation files to annotate microbial network nodes
#'
#' @param ps phyloseq Object, could contains OTU tables, tax table and map table, represented sequences.
#' @param N number object, filter OTU tables by abundance.
#' @param r.threshold Correlation threshold
#' @param p.threshold Significance threshold
#' @param  method method for Correlation calculation,method="pearson" is the default value. The alternatives to be passed to cor are "spearman" and "kendall".
#' @examples
#' nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


nodeadd = function(plotcord =node,otu_table = otu_table,tax_table = tax_table){


  res = merge(plotcord,tax_table,by = "row.names",all = F)
  dim(res)
  head(res)
  row.names(res) = res$Row.names
  res$Row.names = NULL
  plotcord = res

  xx = data.frame(mean  =rowMeans(otu_table))
  head(xx)
  plotcord = merge(plotcord,xx,by = "row.names",all = FALSE)
  head(plotcord)
  # plotcord$Phylum
  row.names(plotcord) = plotcord$Row.names
  plotcord$Row.names = NULL
  head(plotcord)


  return(plotcord)

}


