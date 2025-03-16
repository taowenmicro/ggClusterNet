# library(ggClusterNet)
# library(phyloseq)
# library(tidyverse)
# result = cor_Big_micro(ps = ps,
#                        N = 1000,
#                        # method.scale = "TMM",
#                        r.threshold=0.5,
#                        p.threshold=0.05,
#                        method = "spearman"
# )
#
#
# #--提取相关矩阵
# cor = result[[1]]


cir.maptree2 = function(
    cor = cor,
    n.top = 20,
    r = 1,
    zoom = 0.7

){


  nt = make_igraph(cor)
  dat = node_properties(nt) %>%
    as.data.frame() %>%
    dplyr::arrange(desc(igraph.degree))

  #  需要判断一下多少个点不在我们的num范围之内
  res = node.edge(
    cor = cor,
    select_layout = TRUE,
    clu_method= "cluster_fast_greedy",
    layout_net = "model_maptree2"
  )

  edg = res[[2]]
  # head(edg)
  # head(dat)

  id.r1 = c()
  id.r2 = c()
  i = 1


  for (i in 1:length(dat$igraph.degree)) {
    id = row.names(dat)[i]
    id1 = edg %>% dplyr::filter( OTU_2 %in% c(id)) %>% .$OTU_1 %>% unique()
    id2 = edg %>% dplyr::filter( OTU_1 %in% c(id)) %>% .$OTU_2 %>% unique()
    id0 = c(id1,id2)

    if (length(id0)>n.top) {
      id0 = id0[1:n.top]
    }

    if (i > 1) {
      id0 = setdiff(id0,id.f)
    }

    if (length(id0) == 0) {
      # print(id)
      id.r2 = c(id.r2,id)
    } else{
      id.r1 = c(id.r1,id)
    }

    if (i == 1) {
      id.f = c(id0,id)
    } else{
      id.f = c(id.f,id0,id)
    }
  }


  id.r1 %>% unique() %>% length() +
    id.r2 %>% unique()%>% length()
  unique(c(id.r1,id.r2)) %>% length()


  # [1:num]
  num2 = c(row.names(dat))  %>% unique() %>% length()

  #  中心点和周围点去除之后还有没有点
  num3 = setdiff(row.names(dat),c(unique(id.f,id.r1)))


  dat2 = square_points (num_points = length(id.r1), spacing = 3)
  # ggplot(dat2) + geom_point(aes(x,y)) +theme_void()

  dat2$node = id.r1
  head(dat2)


  for (i in 1:length(id.r1)) {
    id = id.r1[i]
    id1 = edg %>% dplyr::filter( OTU_2 %in% c(id)) %>% .$OTU_1 %>% unique()
    id2 = edg %>% dplyr::filter( OTU_1 %in% c(id)) %>% .$OTU_2 %>% unique()
    id0 = c(id1,id2)

    if (length(id0)>n.top) {
      id0 = id0[1:n.top]
    }

    if (i > 1) {
      id0 = setdiff(id0,id.f)
    }

    if (length(id0) == 0) {

    } else{

      n = length(id0)
      # dat2[i,4]/2
      dat1 = sz.node (x.0 = dat2[i,1],y.0 = dat2[i,2],r = 1,n = n)
      dat1$node = id0
      dat1$class = id
      if (i == 1) {
        dat0 = dat1
        id.f = c(id0,id)
      }else{
        dat0 = rbind(dat0,dat1)
        id.f = c(id.f,id0,id)
      }

    }

  }


  tem = dat0$class %>% table() %>% as.data.frame()
  tem0 = tem %>% dplyr::arrange(desc(Freq))
  head(tem0)
  tem0$. %>% as.character()

  dat2$node = tem0$. %>% as.character()



  for (i in 1:length(id.r1)) {
    id = id.r1[i]
    id1 = edg %>% dplyr::filter( OTU_2 %in% c(id)) %>% .$OTU_1 %>% unique()
    id2 = edg %>% dplyr::filter( OTU_1 %in% c(id)) %>% .$OTU_2 %>% unique()
    id0 = c(id1,id2)
    id.f = c(id.f,id)
    if (length(id0)>n.top) {
      id0 = id0[1:n.top]
    }

    tem = dat2 %>% dplyr::filter(node == id)


    if (i > 1) {
      id0 = setdiff(id0,id.f)
    }

    if (length(id0) == 0) {

    } else{

      n = length(id0)
      # tem[1,4]
      dat1 = sz.node (x.0 = tem[1,1],y.0 = tem[1,2],r =1 ,n = n)
      dat1$node = id0
      dat1$class = id
      if (i == 1) {
        dat0 = dat1
        id.f = c(id0,id)
      }else{
        dat0 = rbind(dat0,dat1)

      }

    }

    id.f = c(id.f,id0)

  }



  intersect(dat0$class,id.r1)
  intersect(dat0$node,id.r1)


  dat0$class2 = "leaf"
  dat2$class2 = "brach"
  dat3 = dat0 %>% full_join(dat2)
  head(dat3)
  colnames(dat3)[1:3] = c("X1","X2","elements")

  head(dat0)
  dat0$class %>% unique()
  dat01 = dat0 %>%dplyr::filter(!node %in% intersect(dat0$node,dat0$class))
  dat02 = dat0 %>%dplyr::filter(node %in% intersect(dat0$node,dat0$class))

  dat0$class%>% unique()%>% length()
  dat02$class %>% unique() %>% length()#
  dat01$class %>% unique() %>% length()#占据大部分，但是少一点





  edge =  data.frame(from = c(dat01$class,dat02$class),to = c(dat01$node,paste0("XXX_" ,1:length(dat02$class))))
  head(edge)
  edge$from %>% unique() %>% length()
  dat2$node %>% length()

  vertices_t  <-  data.frame(
    name = unique(c(as.character(edge$from), as.character(edge$to))))
  head(vertices_t)
  vertices_t$size = sample(1:10,nrow(vertices_t),replace = TRUE)
  intersect(edge$from,edge$to)
  edge$from %>% unique() %>% length()
  mygraph <- igraph::graph_from_data_frame(edge, vertices= vertices_t )

  data = ggraph::create_layout(mygraph, layout = 'circlepack',weight = size)
  head(data)
  # library(ggraph)
  ggraph(mygraph, layout = 'circlepack', sort.by = NULL, direction = "out") +
    geom_node_circle()

  # node = data %>% dplyr::filter(leaf == TRUE ) %>%
  #   dplyr::select(x,y,name)
  # colnames(node) = c("X1","X2", "elements")
  # row.names(node) = node$elements
  # head(data)
  branch = data %>% dplyr::filter(name %in% dat2$node ) %>%
    dplyr::select(x,y,name,r)
  colnames(branch) = c("X1","X2", "elements","r")
  row.names(branch) = branch$elements
  colnames(branch)[1:2] = c("x","y")
  branch$elements = gsub("model_","",branch$elements)
  row.names(branch) = branch$elements
  head(branch)
  branch$elements %>% length()
  #  接入点#-----------
  dat2 = data.frame(x = branch$x,y = branch$y,node = branch$elements,r = branch$r)
  # num5 = dat$igraph.degree[1:num] %>% length()
  dat0$class %>%unique()

  for (i in 1:length(id.r1)) {
    id = id.r1[i]
    id1 = edg %>% dplyr::filter( OTU_2 %in% c(id)) %>% .$OTU_1 %>% unique()
    id2 = edg %>% dplyr::filter( OTU_1 %in% c(id)) %>% .$OTU_2 %>% unique()
    id0 = c(id1,id2)

    if (length(id0)>n.top) {
      id0 = id0[1:n.top]
    }

    if (i > 1) {
      id0 = setdiff(id0,id.f)
    }

    if (length(id0) == 0) {

    } else{

      n = length(id0)
      # dat2[i,4]/2
      dat1 = sz.node (x.0 = dat2[i,1],y.0 = dat2[i,2],r = dat2[i,4]*zoom,n = n)
      dat1$node = id0
      dat1$class = id
      if (i == 1) {
        dat0 = dat1
        id.f = c(id0,id)
      }else{
        dat0 = rbind(dat0,dat1)
        id.f = c(id.f,id0,id)
      }

    }

  }

  head(dat0)




  # ggplot(dat0) + geom_point(aes(x,y))
  # head(dat2)

  #
  # tem = dat0$class %>% table() %>% as.data.frame()
  # tem0 = tem %>% dplyr::arrange(desc(Freq))
  # head(tem0)
  # tem0$. %>% as.character()
  # dat2$node = tem0$. %>% as.character()

  # ggplot(dat2) + geom_point(aes(x,y))


  for (i in 1:length(id.r1)) {
    id = id.r1[i]
    id1 = edg %>% dplyr::filter( OTU_2 %in% c(id)) %>% .$OTU_1 %>% unique()
    id2 = edg %>% dplyr::filter( OTU_1 %in% c(id)) %>% .$OTU_2 %>% unique()
    id0 = c(id1,id2)
    id.f = c(id.f,id)
    if (length(id0)>n.top) {
      id0 = id0[1:n.top]
    }

    tem = dat2 %>% dplyr::filter(node == id)


    if (i > 1) {
      id0 = setdiff(id0,id.f)
    }

    if (length(id0) == 0) {

    } else{

      n = length(id0)
      # tem[1,4]
      dat1 = sz.node (x.0 = tem[1,1],y.0 = tem[1,2],r =tem[1,4]*zoom ,n = n)
      dat1$node = id0
      dat1$class = id
      if (i == 1) {
        dat0 = dat1
        id.f = c(id0,id)
      }else{
        dat0 = rbind(dat0,dat1)

      }

    }

    id.f = c(id.f,id0)

  }



  # intersect(dat0$class,id.r1)
  # intersect(dat0$node,id.r1)


  dat0$class2 = "leaf"
  dat2$class2 = "brach"
  dat3 = dat0 %>% full_join(dat2)
  head(dat3)
  colnames(dat3)[1:3] = c("X1","X2","elements")
  return(dat3)
}


