#' Use the layout in the sne package to calculate the visual layout of the network
#'
#' @param igraph network
#' @examples
#' # Import network data
#' data(igraph)
#' # Computing network related properties：the total number of edges, the number of positive correlation edges, the number of negative correlation edges, the number of vertices...
#' reslt = net_properties(igraph)
#' # export file
#' write.csv(reslt,"netproperties.csv")
#'
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export





net_properties <-function(igraph){
  # igraph.weight <- E(igraph)$weight
  # network property
  # The size of the graph (number of edges)
  num.edges <- length(E(igraph)) # length(curve_multiple(igraph))
  num.edges
  #  Order (number of vertices) of a graph
  num.vertices <- length(V(igraph))# length(diversity(igraph, weights = NULL, vids = 	V(igraph)))
  num.vertices
  #
  connectance <- edge_density(igraph,loops=FALSE)# 同 graph.density;loops如果为TRUE,允许自身环（self loops即A--A或B--B）的存在
  connectance
  # (Average degree)
  average.degree <- mean(igraph::degree(igraph))# 或者为2M/N,其中M 和N 分别表示网络的边数和节点数。
  average.degree
  # (Average path length)
  if (!is.null(E(igraph)$weight)) {
    igraph.weight <- E(igraph)$weight
    E(igraph)$weight = abs(E(igraph)$weight)
  }
  average.path.length <- average.path.length(igraph) # 同mean_distance(igraph) # mean_distance calculates the average path length in a graph
  average.path.length

  # (Diameter)
  diameter <- diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
  diameter


  if (!is.null(E(igraph)$weight)) {
    E(igraph)$weight = igraph.weight
  }
  #  edge connectivity / group adhesion
  edge.connectivity <- edge_connectivity(igraph)
  edge.connectivity
  # (Clustering coefficient)
  clustering.coefficient <- transitivity(igraph)
  clustering.coefficient
  no.clusters <- no.clusters(igraph)
  no.clusters
  # (Degree centralization)
  centralization.degree <- centralization.degree(igraph)$centralization
  centralization.degree
  # (Betweenness centralization)
  centralization.betweenness <- centralization.betweenness(igraph)$centralization
  centralization.betweenness
  # (Closeness centralization)
  centralization.closeness <- centralization.closeness(igraph)$centralization
  centralization.closeness
  if (!is.null(E(igraph)$weight)) {
    num.pos.edges<-sum(igraph.weight>0)# number of postive correlation
    num.neg.edges<-sum(igraph.weight<0)# number of negative correlation
  }else{
    num.pos.edges<-0# number of postive correlation
    num.neg.edges<-0# number of negative correlation
  }


  igraph.network.pro <- rbind(num.edges,num.pos.edges,num.neg.edges,num.vertices,connectance,average.degree,average.path.length,diameter,edge.connectivity,clustering.coefficient,no.clusters,centralization.degree,centralization.betweenness,centralization.closeness)
  rownames(igraph.network.pro)<-c("num.edges","num.pos.edges","num.neg.edges","num.vertices","connectance","average.degree","average.path.length","diameter","edge.connectivity","clustering.coefficient","no.clusters","centralization.degree","centralization.betweenness","centralization.closeness")
  colnames(igraph.network.pro)<- "value"
  return(igraph.network.pro)
}
