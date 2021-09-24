#' Construct edge files, add weights, positive and negative correlation and other basic information
#'
#' @param cor Correlation matrix
#' @param node Node file, containing calculated node coordinates
#' @examples
#' data(ps)
#' result = corMicro (ps = ps,N = 0.02,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # building the node group
#' netClu = data.frame(ID = row.names(cor),group =rep(1:3,length(row.names(cor)))[1:length(row.names(cor))] )
#' netClu$group = as.factor(netClu$group)
#' #Calculate node location
#' result2 = PolygonRrClusterG(cor = cor,nodeGroup = netClu )
#' node = result2[[1]]
#' edge = edgeBuild(cor = cor,plotcord = node)
#'
#' @return edge which contains OTU and its coordinates, the correlation of edge between nodes
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export



edgeBuild = function(cor = cor,plotcord = node){
  # cor <- cor[match( row.names(cor),node$elements),match( row.names(cor),node$elements)]
  dim(cor)

  cor <- cor[match( plotcord$elements,row.names(cor)),match(plotcord$elements, row.names(cor))]
  # colnames(cor)
  #----Use correlation matrix to build network--network package to build network#-----
  g <- network::network(cor, directed=FALSE)
  g

  #---Correlation matrix converted to 0-1
  # origm  <- network::as.matrix.network.adjacency(g)  # get sociomatrix
  edglist <- as.matrix.network.edgelist(g)
  edglist = as.data.frame(edglist)
  #Add weight#---
  set.edge.value(g,"weigt",cor)
  edglist$weight = as.numeric(network::get.edge.attribute(g,"weigt"))
  edges <- data.frame(plotcord[edglist[, 1], ], plotcord[edglist[, 2], ])
  edges$weight = as.numeric(network::get.edge.attribute(g,"weigt"))

  #Here, the edge weights are divided into two categories according to positive and negative
  aaa = rep("a",length(edges$weight))
  for (i in 1:length(edges$weight)) {
    if (edges$weight[i]> 0) {
      aaa[i] = "+"
    }
    if (edges$weight[i]< 0) {
      aaa[i] = "-"
    }
  }
  #add to dges
  edges$wei_label = as.factor(aaa)
  colnames(edges) <- c("X1", "Y1","OTU_1", "X2", "Y2","OTU_2","weight","wei_label")
  edges$midX <- (edges$X1 + edges$X2)/2
  edges$midY <- (edges$Y1 + edges$Y2)/2
 head(edges)
 edges = edges %>% filter(wei_label != "a")

  return(edges)
}

