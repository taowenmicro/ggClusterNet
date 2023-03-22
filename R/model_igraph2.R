#' # The algorithm igraph network layout. while, it saves time and the calculation speed is faster,This is the upgraded version.
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
#' result2 <- model_igraph2(cor = cor,
#' method = "cluster_fast_greedy",
#' seed = 12
#' )
#' node = result2[[1]]
#' @return data.frame
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


model_igraph2 = function(
    cor = cor,
    method = "cluster_fast_greedy",
    seed = 12,
    Top_M = 20

){

  igraph <-  igraph::graph.adjacency(cor, weighted = TRUE, mode = 'undirected')
  # 删除自相关
  igraph <- igraph::simplify(igraph)
  # 删除孤立节点
  igraph <- igraph::delete.vertices(igraph, which(igraph::degree(igraph)==0) )
  # igraph = g1

  # 删除自相关
  igraph <- igraph::simplify(igraph)
  # # 删除孤立节点
  # g <- delete.vertices(g, which(degree(g)==0) )
  col_g <- "#C1C1C1"

  cols <-  colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(Top_M)
  ## 设置网络的weight，为计算模块性做准备
  igraph::E(igraph)$correlation <- igraph::E(igraph)$weight
  igraph::E(igraph)$weight <- abs(igraph::E(igraph)$weight)

  ## 计算网络模块
  set.seed(seed)
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

  igraph::V(igraph)$modularity <- igraph::membership(fc)
  igraph::V(igraph)$label <- igraph::V(igraph)$name
  igraph::V(igraph)$label <- NA
  modu_sort <- igraph::V(igraph)$modularity %>% table() %>% sort(decreasing = T)
  top_num <- Top_M
  modu_name <- names(modu_sort[1:Top_M])
  modu_cols <- cols[1:length(modu_name)]
  names(modu_cols) <- modu_name
  igraph::V(igraph)$color <- igraph::V(igraph)$modularity
  igraph::V(igraph)$color[!(igraph::V(igraph)$color %in% modu_name)] <- col_g
  igraph::V(igraph)$color[(igraph::V(igraph)$color %in% modu_name)] <- modu_cols[match(igraph::V(igraph)$color[(igraph::V(igraph)$color %in% modu_name)],modu_name)]
  igraph::V(igraph)$frame.color <- igraph::V(igraph)$color

  E(igraph)$color <- col_g
  for ( i in modu_name){
    col_edge <- cols[which(modu_name==i)]
    otu_same_modu <-igraph::V(igraph)$name[which(igraph::V(igraph)$modularity==i)]
    igraph::E(igraph)$color[(data.frame(igraph::as_edgelist(igraph))$X1 %in% otu_same_modu)&(data.frame(igraph::as_edgelist(igraph))$X2 %in% otu_same_modu)] <- col_edge
  }



  # 计算 layout
  sub_net_layout <- igraph::layout_with_fr(igraph, niter=999,grid = 'nogrid')
  data = as.data.frame(sub_net_layout)
  data$OTU = igraph::get.vertex.attribute(igraph)$name
  colnames(data) = c("X1","X2","elements")

  tem =  igraph::V(igraph)$modularity
  tem[!tem %in% modu_name] = "mini_model"
  tem[tem %in% modu_name] = paste("model_",tem[tem %in% modu_name],sep = "")

  row.names(data) = data$elements
  dat = data.frame(orig_model = igraph::V(igraph)$modularity,
                   model = tem,
                   color = igraph::V(igraph)$color,
                   OTU = igraph::get.vertex.attribute(igraph)$name,
                   X1 = data$X1,
                   X2 = data$X2
  )

  return(list(data,dat,igraph))


}


