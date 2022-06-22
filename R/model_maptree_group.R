

model_maptree_group = function(
    cor = cor,
    nodeGroup =netClu,
    seed = 12
  ){
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
  # Output igraph object
  igraph  = igraph::graph_from_data_frame(edges, directed = FALSE, vertices = node)

  #-Extraction module
  # nodeGroup = data.frame(ID = names(membership(fc)),group = as.vector(membership(fc)))
  # table(nodeGroup$group)
  result4 = nodeEdge(cor = cor)
  edge = result4[[1]]
  node = result4[[2]]
  igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
  igraph.degree<-igraph::degree(igraph) %>% as.data.frame()
  colnames(igraph.degree) = "degree"
  igraph.degree$ID = row.names(igraph.degree)
  nodeGroup <- nodeGroup %>% left_join(igraph.degree,na_matches = "never")
  nodeGroup$degree[is.na(nodeGroup$degree)] = 0
  nodeGroup <- nodeGroup %>%
    dplyr::arrange(desc(degree))
  head(nodeGroup)

  edge =  data.frame(model = paste("model_",nodeGroup$group,sep = ""),OTU = nodeGroup$ID)
  head(edge)
  colnames(edge) = c("from","to")

  vertices_t  <-  data.frame(
    name = unique(c(as.character(edge$from), as.character(edge$to))))
  head(vertices_t)
  vertices_t$size = sample(1:10,nrow(vertices_t),replace = TRUE)

  mygraph <- igraph::graph_from_data_frame(edge, vertices= vertices_t )
  #-----------------------------------设置颜色映射参数-------------------------
  # ?create_layout
  # ,weight = mean, sort.by = NULL, direction = "out"
  data = ggraph::create_layout(mygraph, layout = 'circlepack',weight = size)
  head(data)

  node = data %>% dplyr::filter(leaf == TRUE ) %>%
    dplyr::select(x,y,name)
  colnames(node) = c("X1","X2", "elements")
  row.names(node) = node$elements
  branch = data %>% dplyr::filter(leaf != TRUE ) %>%
    dplyr::select(x,y,name)
  colnames(branch) = c("X1","X2", "elements")
  row.names(branch) = branch$elements
  colnames(branch)[1:2] = c("x","y")
  branch$elements = gsub("model_","",branch$elements)
  row.names(branch) = branch$elements

  return(list(node,branch))
}

