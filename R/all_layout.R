

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
    node = total_layout(layout =layout,method = "cluster_fast_greedy" )

  } else{

    #--make groups
    if (model == TURE) {
      netClu  = modulGroup( cor = cor,cut = 3,method = method )
    } else if(model != TURE){
      netClu = group
    }

    #group layout
    node = group_layout(layout =layout,
                            netClu = netClu,
                            sna_layout = "circle"
    )



  }

  return(node)

}

total_layout = function(layout =layout,method = "cluster_fast_greedy" ){

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
  #3
  if (layout == "model_maptree") {
    result2 <- model_maptree(cor = cor,
                             method = method,
                             seed = 12
    )
    node = result2
  }

 return(node)


}


group_layout = function(layout =layout,
                        netClu = netClu,
                        sna_layout = "circle"
                        ){

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
  return(node)

}
