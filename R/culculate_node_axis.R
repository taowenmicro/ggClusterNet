


#' Call the network layout algorithm for coordinate calculation
#'
#' @title use layout calculated the axis a and y
#' @description filter microbiome data
#' @param cor.matrix Correlation matrix
#' @param  layout network layout algorithm in ggClusterNet
#' @examples
#' data(cor)
#'node = culculate_node_axis(
#'  cor.matrix = cor,
#'  layout = "model_Gephi.2",
#'  seed = 1,
#'  group = NULL,
#'  model = FALSE,
#'  method = "cluster_fast_greedy")

culculate_node_axis = function(
  cor.matrix = cor,
  layout = "model_Gephi.2",
  seed = 1,
  group = NULL,
  model = FALSE,
  method = "cluster_fast_greedy"
){

  if (is.null(group)) {
    #--total layout
    node = total_layout(cor = cor.matrix,layout =layout,method = method )
  }else{

    if (model) {

      netClu  = modulGroup( cor = cor.matrix,cut = 3,method = method )

      }else{
      netClu = group
      }


    node = group_layout(
      cor = cor.matrix,
      layout =layout,
                            netClu = netClu,
                            sna_layout = "circle"
    )

  }
  return(node)
}

total_layout = function(cor = cor,layout =layout,method = "cluster_fast_greedy" ){

  #--1
  if (layout == "model_Gephi.2") {
    result2 <- model_Gephi.2(cor = cor,
                             method = method,
                             seed = 12
    )
    node = result2[[1]]
  }
  #2
  if (layout == "model_igraph") {
    result2 <- model_igraph(cor = cor,
                            method = method,
                            seed = 12
    )
    node = result2[[1]]
  }

  if (layout == "model_igraph2") {
    result2 <- model_igraph2(cor = cor,
                            method = method,
                            seed = 12
    )
    node = result2[[1]]
  }

  #3
  if (layout == "model_maptree") {
    result2 <- model_maptree(cor = cor,
                             method = method,
                             seed = 12
    )
    node = result2[[1]]
  }

  if (layout == "model_maptree2") {
    result2 <- model_maptree2(cor = cor,
                             method = method,
                             seed = 12
    )
    node = result2[[1]]
  }

  if (layout == "model_Gephi.3") {
    result2 = model_Gephi.3(
      cor = cor,
      t0 = 0.98,
      t2 = 1.2,
      t3 = 0)

    node0 = result2[[1]]
    node=data.frame(row.names = node0$OTU,X1 = node0$X1,X2 = node0$X2,elements = node0$OTU )

  }

  if (layout == "model_Gephi.4") {

  dat <- layout.forceatlas2(cor = cor,
                            noigraph = TRUE,
                            iterations=500,
                            plotstep=10,
                            directed=FALSE,
                            linlog = F,
                            pos = NULL,
                            nohubs = F,
                            k = 100,
                            gravity=0.1,
                            ks=0.1,
                            ksmax=10,
                            delta = 15,
                            center=NULL,
                            tolerance = 5,
                            dim = 2,
                            plotlabels=TRUE
  )



  # head(dat)
  node =dat
  row.names(node) = node$elements
  }




 return(node)


}


group_layout = function(
  cor = cor,
  layout =layout,
  netClu = netClu,
  sna_layout = "circle"){
  if (layout == "randomClusterG") {
    result2 = randomClusterG (cor = cor,nodeGroup =netClu )
    node = result2[[1]]
    head(node)
  }

  if (layout == "PolygonClusterG") {
    result2 = PolygonClusterG(cor = cor,nodeGroup =netClu )
    node = result2[[1]]
    head(node)
  }

  if (layout == "PolygonRrClusterG") {
    #
    result2 = PolygonRrClusterG(cor = cor,nodeGroup =netClu,zoom = 1,zoom2 = 1)
    node = result2[[1]]
    head(node)
  }

  if (layout == "ranSNEClusterG") {
    result2 = ranSNEClusterG (cor=  cor,layout = sna_layout)
    node = result2[[1]]
    head(node)
  }

  if (layout == "PolygonModsquareG") {
    result2 <- PolygonModsquareG(cor = cor,nodeGroup =netClu,r1 = 1,N = 0.5,cut = 5)
    node = result2[[1]]
  }

  if (layout == "PolyRdmNotdCirG") {
    result2 = PolyRdmNotdCirG (cor = cor,nodeGroup =netClu )
    node = result2[[1]]
    head(node)
  }
  if (layout == "model_filled_circle") {
    #-实心圆2
    result2 = model_filled_circle(cor = cor,
                                  culxy =TRUE,
                                  da = NULL,# 数据框，包含x,和y列
                                  nodeGroup =netClu,
                                  mi.size = 1,# 最小圆圈的半径，越大半径越大
                                  zoom = 0.15# 不同模块之间距离
                                  )
    node = result2[[1]]
    head(node)
  }

  if (layout == "model_maptree_group") {
    result2 = model_maptree_group (
      cor = cor,
      nodeGroup =netClu,
      seed = 12)
    node = result2[[1]]
    head(node)
  }

  return(node)

}
