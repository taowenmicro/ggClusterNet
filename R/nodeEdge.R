
#' Construct igraph and tidygraph object
#'
#' @param cor Correlation matrix
#' @param zero  TRUE or FALSE. If the edge-to-edge correlation is 0, FALSE preserves the output.
#' @examples
#' data(ps)
#' result = corMicro (ps = ps,N = 0.02,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # nodeEdge function generates nodes and edge files
#' result4 = nodeEdge(cor = cor)
#' # Extract edge file
#' edge = result4[[1]]
#' dim(edge)
#' edge$weight
#' # Extract node file
#' node = result4[[2]]
#' dim(node)
#' # Input to the igraph for plot
#' igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
#' # Input to the tidygraph for plot
#' library(tidygraph)
#' library(ggraph)
#' tbl_graph = tidygraph::tbl_graph(nodes = node, edges = edge, directed = FALSE)
#' tbl_graph
#' @return list Documents required for drawing
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export

nodeEdge = function(corr = cor,zero = TRUE){

  # 构造边文件
  edges <- tibble::tibble(from = rep(row.names(corr), ncol(corr)),
                          to = rep(colnames(corr), each = nrow(corr)),
                          r = as.vector(corr)
  )
  head(edges)
  # 提取一半的矩阵，并去除对角线的相关（自己于自己相关）
  edges <- dplyr::filter(edges, as.vector(lower.tri(corr)))
  #按照R 阈值和P阈值筛选
if (zero == TRUE) {
  edges <- dplyr::filter(edges,abs(r)> 0)
  edges$r
  edgesn = edges
}

  if (zero == FALSE) {
    edgesn = edges
    edges <- dplyr::filter(edges,abs(r)> 0)

  }
  colnames(edges)[3] = "weight"
  #---设置边的正负
  # E.color <- edges$weight
  edges$direction <- ifelse(edges$weight>0, "pp",ifelse(edges$weight<0, "np","ns"))
  dim(edges)

  # nodes<- dplyr::select(nodes, elements, everything())
  # colnames(edges)[1] = c("name")

  # unique(nodes$elements)

  node = tibble::tibble(name = unique(c(edgesn$from, edgesn$to)))
  row.names(node) = node$name
  dim(node)




  return(list(edges,node))
}









