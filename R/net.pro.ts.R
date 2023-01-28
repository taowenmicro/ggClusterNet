
net.pro.ts = function(
    dat = dat
){
  for (i in 1:length(names(dat))) {
    cor = dat[[names(dat)[i]]]
    igraph  = igraph::graph_from_data_frame(nodeEdge(cor = cor)[[1]],
                                            directed = FALSE,
                                            vertices = nodeEdge(cor = cor)[[2]])
    # nodepro = node_properties(igraph)
    netpro_result<- net_properties.4(igraph)
    colnames(netpro_result) = names(dat)[i]

    if (i == 1 ) {
      tab = netpro_result
    } else{
      tab = cbind(tab,netpro_result)
    }
  }
  return(tab)
}
