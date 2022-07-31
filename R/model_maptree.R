#' # The algorithm imitate maptree network layout. while, it saves time and the calculation speed is faster,This is the upgraded version.
#'
#' @title This is the upgraded version.model_Gephi.2:The algorithm imitate Gephi's network layout. while, it saves time and the calculation speed is faster
#' @description Enter correlation matrix, calculate network modules, and Calculate the coordinates of the node.
#' @param cor Correlation matrix Clustering Algorithm
#' @param method Clustering Algorithm
#' @param seed Set random seed
#' @details
#' By default, returns a list
#' \itemize{
#' \item{cluster_fast_greedy: }
#' \item{cluster_walktrap: }
#' \item{cluster_edge_betweenness: }
#' \item{cluster_spinglass: }
#' }
#' @examples
#' data(ps)
#' result = corMicro (ps = ps,N = 100,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # building the node group
#' result2 <- model_maptree(cor = cor,
#' method = "cluster_fast_greedy",
#' seed = 12
#' )
#' node = result2
#' @return data.frame
#' @author Contact: Tao Wen \email{taowen@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn} yongxin liu \email{yxliu@@genetics.ac.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Tao Wen#, Penghao Xie#, Shengdie Yang, Guoqing Niu, Xiaoyu Liu, Zhexu Ding, Chao Xue, Yong-Xin Liu *, Qirong Shen, Jun Yuan*
#' ggClusterNet: an R package for microbiome network analysis and modularity-based multiple network layouts
#' iMeta 2022,DOI: \url{doi: 10.1002/imt2.32}
#' @export


model_maptree = function(
  cor = cor,
  method = "cluster_fast_greedy",
  seed = 12
){
  corr <- cor

  # Construct Edge File
  edges <- data.frame(from = rep(row.names(corr), ncol(corr)),
                      to = rep(colnames(corr), each = nrow(corr)),
                      r = as.vector(corr)
  )
  # Extract half of the matrix, and remove the diagonal correlation (self related to yourself)
  edges <- dplyr::filter(edges, as.vector(lower.tri(corr)))
  colnames(edges)[3] = "weight"
  #---Set the sign of the edge
  # E.color <- edges$weight
  edges$direction <- ifelse(edges$weight>0, "pp",ifelse(edges$weight<0, "np","ns"))
  node = data.frame(name = unique(c(as.character(edges$from),as.character( edges$to))))
  row.names(node) = node$name
  # Output igraph object
  # igraph  = make_igraph(corr)
  igraph  = igraph::graph_from_data_frame(edges, directed = FALSE, vertices = node)
  igraph::E(igraph)
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

  modularity <- modularity(igraph,membership(fc))
  #-Extraction module
  netClu = data.frame(ID = names(membership(fc)),group = as.vector(membership(fc)))
  table(netClu$group)
  result4 = nodeEdge(cor = cor)
  edge = result4[[1]]
  node = result4[[2]]
  igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
  igraph.degree<-igraph::degree(igraph) %>% as.data.frame()
  colnames(igraph.degree) = "degree"
  igraph.degree$ID = row.names(igraph.degree)
  netClu <- netClu %>% left_join(igraph.degree,na_matches = "never")
  netClu$degree[is.na(netClu$degree)] = 0
  netClu <- netClu %>%
    dplyr::arrange(desc(degree))
  head(netClu)

  edge =  data.frame(model = paste("model_",netClu$group,sep = ""),OTU = netClu$ID)
  head(edge)
  colnames(edge) = c("from","to")

  vertices_t  <-  data.frame(
    name = unique(c(as.character(edge$from), as.character(edge$to))))
  head(vertices_t)
  vertices_t$size = sample(1:10,nrow(vertices_t),replace = TRUE)

  mygraph <- igraph::graph_from_data_frame(edge, vertices= vertices_t )
  #-----------------------------------设置颜色映射参数-------------------------
  # ?create_layout
  # ,weight = mean, sort.by = NULL, direction = "out"
  data = ggraph::create_layout(mygraph, layout = 'circlepack',weight = size)
  head(data)

  node = data %>% dplyr::filter(leaf == TRUE ) %>%
    dplyr::select(x,y,name)
  colnames(node) = c("X1","X2", "elements")
  row.names(node) = node$elements
  return(list(node,netClu))
}

