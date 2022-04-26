#' # This is the function named 'modulGroup'# which Enter the correlation matrix. There are four module clustering algorithms inside, which can calculate the modularity and then group according to the modularity.
#'
#' @title Correlation matrix calculates groups according to the degree of network module speech.
#' @description Enter correlation matrix, calculate network modules, and generate groups.
#' @param cor Correlation matrix
#' @param nodeGroup Microbiome grouping table
#' @param seed set random seed,Default 12
#' @details
#' By default, returns table, contain node and group imformation
#' The available method to culculate Degree of modularity include the following:
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
#' PolygonMaptreeG(cor = cor,nodeGroup = gru,seed = 12)
#'
#' @return data table.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


PolygonMaptreeG = function(
  cor = cor,
  nodeGroup = gru,
  seed = 12
){
  diag(cor) = 0
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

  #-Extraction module
  netClu =  nodeGroup
  table(netClu$group)
  result4 = nodeEdge(cor = cor)
  edge = result4[[1]]
  node = result4[[2]]
  igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
  igraph.degree<-igraph::degree(igraph) %>% as.data.frame()
  colnames(igraph.degree) = "degree"
  igraph.degree$ID = row.names(igraph.degree)
  netClu <- netClu %>%
    dplyr::left_join(igraph.degree,na_matches = "never")
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

  data = ggraph::create_layout(mygraph, layout = 'circlepack',weight = size)
  head(data)

  node = data %>% dplyr::filter(leaf == TRUE ) %>%
    dplyr::select(x,y,name)
  colnames(node) = c("X1","X2", "elements")
  row.names(node) = node$elements
  return(node)
}

