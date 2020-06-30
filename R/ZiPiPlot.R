#' Use the layout in the sne package to calculate the visual layout of the network
#'
#' @param cor Correlation matrix
#' @param method method to culculate Degree of modularity
#' @examples
#' data(igraph)
#' res = ZiPiPlot(igraph = igraph,method = "cluster_fast_greedy")
#' res[[1]]
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


ZiPiPlot = function(igraph = igraph,method = "cluster_fast_greedy"){




  if (method == "cluster_walktrap" ) {
    fc <- igraph::cluster_walktrap(igraph,weights =  abs(E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }

  if (method == "cluster_edge_betweenness" ) {
    fc <- igraph::cluster_edge_betweenness(igraph,weights =  abs(E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }
  if (method == "cluster_fast_greedy" ) {
    fc <- igraph::cluster_fast_greedy(igraph,weights =  abs(E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }
  if (method == "cluster_spinglass" ) {
    fc <- igraph::cluster_spinglass(igraph,weights =  abs(E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }
  modularity <- igraph::modularity(igraph,membership(fc))
  # 模块化程度
  modularity
  # 按照模块为节点配色，这里我们可以加入到nodes中
  comps <- membership(fc)

  # comps
  V(igraph)$module <- as.character(comps)

  taxa.roles <- microbiomeSeq::module.roles(igraph)

  head(taxa.roles)
  head(taxa.roles)
  taxa.roles$label = row.names(taxa.roles)
  for (i in 1:nrow(taxa.roles))if(taxa.roles[i,3]> 0.62|taxa.roles[i,1]> 2.5) {
    taxa.roles[i,5]=taxa.roles[i,5]
  }else{
    taxa.roles[i,5]= ""
  }

  library(ggrepel)
  p <- microbiomeSeq::plot_roles(taxa.roles) + ggrepel::geom_text_repel(data = taxa.roles, aes(x = p, y = z, color = module,label=taxa.roles$label),size=4)#
  #geom_text(data = taxa.roles, aes(x = p, y = z, color = module,label=taxa.roles$label),size=4)
  # print(p)

  return(list(p,taxa.roles))
}

