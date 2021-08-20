#' # This is the function named 'modulGroup'# which Enter the correlation matrix. There are four module clustering algorithms inside, which can calculate the modularity and then group according to the modularity.
#'
#' @title Correlation matrix calculates groups according to the degree of network module speech.
#' @description Enter correlation matrix, calculate network modules, and generate groups.
#' @param cor Correlation matrix
#' @param cut which model contain node less than cut,will be tegether to other group,if NULL was assignment, will not to do with nothing.
#' @param method method to culculate Degree of modularity.There are four module clustering algorithms inside.
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
#' result = corMicro (ps = ps,N = 0.02,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # building the node group
#' netClu  = modulGroup( cor = cor,cut = 3,method = "cluster_fast_greedy" )
#' netClu$group = as.factor(netClu$group)
#' result2 = ranSNEClusterG (cor=  cor,layout ="circle")
#' node = result2[[1]]
#'
#' @return data table.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


modulGroup = function( corr = cor,cut = NULL,method = "cluster_walktrap"){

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


  # nodes<- dplyr::select(nodes, elements, everything())
  # colnames(edges)[1] = c("name")

  # unique(nodes$elements)

  node = data.frame(name = unique(c(as.character(edges$from),as.character( edges$to))))
  row.names(node) = node$name

  # Output igraph object
  igraph  = igraph::graph_from_data_frame(edges, directed = FALSE, vertices = node)

  # Modularity - there are four common calculation methods here


  if (method == "cluster_walktrap" ) {
    fc <- cluster_walktrap(igraph,weights =  abs(E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }

  if (method == "cluster_edge_betweenness" ) {
    fc <- cluster_edge_betweenness(igraph,weights =  abs(E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }
  if (method == "cluster_fast_greedy" ) {
    fc <- cluster_fast_greedy(igraph,weights =  abs(E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }
  if (method == "cluster_spinglass" ) {
    fc <- cluster_spinglass(igraph,weights =  abs(E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }

  modularity <- modularity(igraph,membership(fc))
  # Modularity
  modularity
  #-Extraction module

  netClu = data.frame(ID = names(membership(fc)),group = as.vector(membership(fc)))
  head(netClu)

  #-We see that many modules have only one node, so these cannot be used as a module, we merge and then define as other
  #


  mod = as.data.frame(table(netClu$group))

  # mod$Freq[mod$Freq< cut] = 0
  # All modules with less than one cut of otu are grouped into one, defined as 0

  if (is.null(cut)) {
    netClu = netClu
  } else {
    for (i in 1:length(netClu$group)) {
      if (netClu$group[i] %in% as.character(mod$Var1[mod$Freq>= cut])) {
        netClu$group[i] = netClu$group[i]
      } else {
        netClu$group[i] = 0
      }
    }
  }

  netClu$group = as.factor(netClu$group)
  return(netClu)
}

