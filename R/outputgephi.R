#' use the phylsoeq object and cor matrix building Gephi input network file contain edge and node
#'
#' @param ps phyloseq abject
#' @param cor Correlation matrix
#' @examples
#' data(ps)
#' result = corMicro (ps = ps,N = 0.02,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' gephi <- outputgephi(ps = result[[3]],cor =  cor)
#'edge_Gephi <- gephi[[1]]
#'head(edge_Gephi)
#'node_Gephi <- gephi[[2]]
#'head(node_Gephi)
#' @return result2 Which contains 2 lists.Result2[[1]], consists of OTU and its corresponding coordinates.
#' Result2[[2]], consists of the network center coordinates of each group
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export

outputgephi <- function(ps = ps,cor =  cor){
  #--imput cor matrix
  result4 = nodeEdge(cor = cor)
  edge = result4[[1]]
  node = result4[[2]]
  # edge
  edge_Gephi = data.frame(source = edge$from,target = edge$to,correlation =  edge$weight,direct= "undirected",cor =  edge$direction)
  tax_table = as.data.frame(vegan_tax(ps))
  otu_table = as.data.frame(t(vegan_otu(ps)))
  #---node add tax and mean abundance of all of sample #-----------
  node_Gephi = nodeadd(plotcord =node,otu_table = otu_table ,tax_table = tax_table)
  return(list(edge_Gephi,node_Gephi))
}


