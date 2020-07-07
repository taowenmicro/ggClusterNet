#' Compute network node attributes
#'
#' @param igraph network
#' @examples
#' data(igraph)
#' # Computing node related properties: the number of connections to the node
#' node_properties(igraph)
#' dim(nodepro)
#' write.csv(nodepro,"nodeproperties.csv")
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


node_properties<-function(igraph,outdir){

  # Nodal degree
  igraph.degree<-igraph::degree(igraph)
  # Nodal centrality
  igraph.cen.degree<-centralization.degree(igraph)$res
  #
  igraph.betweenness<-centralization.betweenness(igraph)$res
  #
  igraph.closeness<-centralization.closeness(igraph)$res

  igraph.node.pro <- cbind(igraph.degree,igraph.closeness,igraph.betweenness,igraph.cen.degree)
  colnames(igraph.node.pro)<-c("igraph.degree","igraph.closeness","igraph.betweenness","igraph.cen.degree")
  igraph.node.pro
}


