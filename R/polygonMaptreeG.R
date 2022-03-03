

# PolygonMaptreeG(cor = cor,group = gru,seed = 12)

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

