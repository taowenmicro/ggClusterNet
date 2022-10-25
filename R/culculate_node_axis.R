
# node = culculate_node_axis(
#   cor.matrix = cor,
#   layout = "model_Gephi.2",
#   seed = 1,
#   group = NULL,
#   model = FALSE,
#   method = "cluster_fast_greedy")


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
    node = total_layout(cor = cor.matrix,layout =layout,method = "cluster_fast_greedy" )
  }else{ if (model) {
      netClu  = modulGroup( cor = cor.matrix,cut = 3,method = method )}else{
      netClu = group}
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
    result2 = PolygonRrClusterG(cor = cor,nodeGroup =netClu,zoom1 = 0.5,zoom2 = 0.2)
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
                                  zoom = 0.3# 不同模块之间距离
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
