#' Add annotation information to network nodes
#'
#' @param plotcord The coordinate information of node OTU
#' @param otu_table Relative abundance of node OTU
#' @param tax_table Species annotation information of the node OTU
#' @examples
#' data(ps)
#' # Calculate the correlation network of microbial community data.Generate 4 lists
#' result = corMicro (ps = ps,N = 100,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#'# Extract phyloseq object from result
#' ps_net = result[[3]]
#' otu_table = as.data.frame(t(vegan_otu(ps_net)))
#' Extract tax_table from ps_net
#' vegan_tax <-  function(physeq){
#' tax <-  tax_table(physeq)
#'
#' return(as(tax,"matrix"))
#' }
#' tax_table = as.data.frame(vegan_tax(ps_net))
#' # Set up OTU grouping
#' netClu = data.frame(ID = row.names(otu_table),group =rep(1:5,length(row.names(otu_table)))[1:length(row.names(otu_table))] )
#' netClu$group = as.factor(netClu$group)
#' # Construct a network layout. Arrange network nodes to different locations according to grouping
#' result2 = PolygonClusterG (cor = cor,nodeGroup =netClu  )
#' node = result2[[1]]
#' nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export

nodeadd = function(plotcord =node,otu_table = otu_table,tax_table = tax_table){


  res = base::merge(plotcord,tax_table,by = "row.names",all.x = TRUE)
  dim(res)
  head(res)
  row.names(res) = res$Row.names
  res$Row.names = NULL
  plotcord = res

  xx = data.frame(mean  =rowMeans(otu_table))
  head(xx)
  plotcord =  base::merge(plotcord,xx,by = "row.names",all.x = TRUE)
  head(plotcord)
  # plotcord$Phylum
  row.names(plotcord) = plotcord$Row.names
  plotcord$Row.names = NULL
  head(plotcord)


  return(plotcord)

}


