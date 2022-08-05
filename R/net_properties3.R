
net_properties.3 <-function(igraph, n.hub = FALSE
){
  # igraph.weight <- E(igraph)$weight
  # network property
  # The size of the graph (number of edges)
  num.edges <- length(igraph::E(igraph))
  num.edges
  #  Order (number of vertices) of a graph
  num.vertices <- length(igraph::V(igraph))# length(diversity(igraph, weights = NULL, vids = 	V(igraph)))
  num.vertices
  
  connectance <- igraph::edge_density(igraph,loops=FALSE)
  
  # (Average degree)
  average.degree <- mean(igraph::degree(igraph))# 或者为2M/N,其中M 和N 分别表示网络的边数和节点数。
  average.degree
  # (Average path length)
  if (!is.null(igraph::E(igraph)$weight)) {
    igraph.weight <- igraph::E(igraph)$weight
    igraph::E(igraph)$weight = abs(igraph::E(igraph)$weight)
  }
  
  average.path.length <- igraph::average.path.length(igraph) # 同mean_distance(igraph) # mean_distance calculates the average path length in a graph
  average.path.length
  
  # (Diameter)
  diameter <- igraph::diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
  diameter
  
  
  if (!is.null(igraph::E(igraph)$weight)) {
    igraph::E(igraph)$weight = igraph.weight
  }
  #  edge connectivity / group adhesion
  edge.connectivity <- igraph::edge_connectivity(igraph)
  edge.connectivity
  # (Clustering coefficient)
  clustering.coefficient <- igraph::transitivity(igraph,type = "average")
  clustering.coefficient
  
  no.clusters <- igraph::no.clusters(igraph)
  no.clusters
  # (Degree centralization)
  centralization.degree <- igraph::centralization.degree(igraph)$centralization
  centralization.degree
  # (Betweenness centralization)
  centralization.betweenness <- igraph::centralization.betweenness(igraph)$centralization
  centralization.betweenness
  # (Closeness centralization)
  centralization.closeness <- igraph::centralization.closeness(igraph)$centralization
  centralization.closeness
  if (!is.null(E(igraph)$weight)) {
    num.pos.edges<-sum(igraph.weight>0)# number of postive correlation
    num.neg.edges<-sum(igraph.weight<0)# number of negative correlation
  }else{
    num.pos.edges<-0# number of postive correlation
    num.neg.edges<-0# number of negative correlation
  }
  
  #-----add RM #------
  modularity_igraph = function(net,method = "cluster_walktrap"){
    if (method == "cluster_walktrap" ) {
      fc <- cluster_walktrap(igraph,weights =  abs(igraph::E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
    }
    
    if (method == "cluster_edge_betweenness" ) {
      fc <- cluster_edge_betweenness(igraph,weights =  abs(igraph::E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
    }
    if (method == "cluster_fast_greedy" ) {
      fc <- cluster_fast_greedy(igraph,weights =  abs(igraph::E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
    }
    if (method == "cluster_spinglass" ) {
      fc <- cluster_spinglass(igraph,weights =  abs(igraph::E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
    }
    modularity <- igraph::modularity(net,membership(fc))
    return(modularity)
  }
  
  mod1 = modularity_igraph(igraph,method = "cluster_walktrap")
  rand.g <- erdos.renyi.game(length(V(igraph)), length(E(igraph)),type = "gnm")
  mod2 = modularity_igraph(rand.g,method = "cluster_walktrap")
  RM = (mod1-mod2)/mod2
  #---
  if (n.hub) {
    res = ZiPiPlot(igraph = igraph,method = "cluster_walktrap")
    data = res[[2]]
    head(data)
    n.hub = data$roles[data$roles != "Peripherals" ] %>% length()
  } else {
    n.hub = "Not.calculated"
  }
  
  
  igraph.network.pro <- rbind(num.edges,num.pos.edges,
                              num.neg.edges,num.vertices,
                              connectance,average.degree,average.path.length,
                              diameter,edge.connectivity,clustering.coefficient,
                              no.clusters,
                              centralization.degree,centralization.betweenness,
                              centralization.closeness,
                              RM,
                              n.hub
                              
  )
  rownames(igraph.network.pro)<-c("num.edges(L)",
                                  "num.pos.edges",
                                  "num.neg.edges",
                                  "num.vertices(n)",
                                  "Connectance(edge_density)","average.degree(Average K)","average.path.length",
                                  "diameter","edge.connectivity",
                                  "mean.clustering.coefficient(Average.CC)",
                                  "no.clusters","centralization.degree",
                                  "centralization.betweenness","centralization.closeness",
                                  "RM(relative.modularity)",
                                  "the.number.of.keystone.nodes"
  )
  colnames(igraph.network.pro)<- "value"
  return(igraph.network.pro)
}

# library(ggClusterNet)
# data("igraph")
# library(igraph)
# #-运行一下net_properties.3函数
# # 在进行下面代码运算
# net_properties.3(igraph,T)
