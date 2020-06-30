#' Construct a network layout. Arrange network nodes to different locations according to grouping
#'
#' @param cor Correlation matrix
#' @param zero  TRUE or FEASE
#' @examples
#' data(ps)
#' result = corMicro (ps = ps,N = 0.02,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # building the node group
#' netClu = data.frame(ID = row.names(cor),group =rep(1:3,length(row.names(cor)))[1:length(row.names(cor))] )
#' netClu$group = as.factor(netClu$group)
#' result2 = PolygonRrClusterG (cor = cor,nodeGroup =netClu ）
#' node = result2[[1]]
#'
#'
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


nodeEdge = function(cor = cor,zero = TRUE){
  corr <-cor
  # 构造边文件
  edges <- tibble::tibble(from = rep(row.names(corr), ncol(corr)),
                          to = rep(colnames(corr), each = nrow(corr)),
                          r = as.vector(corr)
  )
  head(edges)
  # 提取一半的矩阵，并去除对角线的相关（自己于自己相关）
  edges <- dplyr::filter(edges, lower.tri(corr))
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
  library(tidyverse)
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












