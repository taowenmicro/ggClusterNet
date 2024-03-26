
node.edge = function(
  cor = cor,
  select_layout = TRUE,
    clu_method=clu_method,
   layout_net = layout_net,
  group.node = NULL,
  model.node = FALSE
){

  if (select_layout) {
    node = NULL
    node = culculate_node_axis(
      cor.matrix = cor,
      layout = layout_net,
      seed = 1,
      group = group.node,
      model = model.node,
      method = clu_method)

  }else {
    result2 <- model_Gephi.2(cor = cor,
                             method = clu_method,
                             seed = 12
    )
    node = result2[[1]]
  }

  head(node)
  edge = edgeBuild(cor = cor,node = node)
  head(edge)
  igraph  = igraph::graph_from_data_frame(nodeEdge(cor = cor)[[1]],
                                          directed = FALSE,
                                          vertices = nodeEdge(cor = cor)[[2]])
  nodepro = node_properties(igraph)
  nodeG = merge(node,nodepro,by = "row.names",all.x  = TRUE)
  row.names(nodeG) = nodeG$Row.names
  nodeG$Row.names = NULL

  numna = (dim(nodeG)[2] - 3) : dim(nodeG)[2]
  nodeG[,numna][is.na(nodeG[,numna])] = 0
  return(list(node = nodeG,edge = edge))
}

