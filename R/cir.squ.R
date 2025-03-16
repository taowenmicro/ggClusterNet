


cir.squ = function(
    cor = cor,
    n.top = 20,
    r = 1

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

  # num5 = dat$igraph.degree[1:num] %>% length()



  for (i in 1:length(id.r1)) {
    id = id.r1[i]
    id1 = edg %>% dplyr::filter( OTU_2 %in% c(id)) %>% .$OTU_1 %>% unique()
    id2 = edg %>% dplyr::filter( OTU_1 %in% c(id)) %>% .$OTU_2 %>% unique()
    id0 = c(id1,id2) %>% unique()

    id0 = setdiff(id0,id.r1)





    if (length(id0)>n.top) {
      id0 = id0[1:c(n.top)]

    }

    if (i > 1) {
      id0 = setdiff(id0,id.f)
    }

    if (length(id0) == 0) {

    } else{

      n = length(id0)
      dat1 = sz.node (x.0 = dat2[i,1],y.0 = dat2[i,2],r = r,n = n)
      dat1$node = id0
      dat1$class = id
      if (i == 1) {
        dat0 = dat1
        id.f = c(id0)
      }else{
        dat0 = rbind(dat0,dat1)
        id.f = c(id.f,id0,id)
      }

    }

  }


  head(dat0)
#
# intersect(dat0$node,dat0$class)
#
#   head(dat2)


  tem = dat0$class %>% table() %>% as.data.frame()
  tem0 = tem %>% dplyr::arrange(desc(Freq))
  head(tem0)
  tem0$. %>% as.character()

  dat2$node = tem0$. %>% as.character()

  # ggplot(dat2) + geom_point(aes(x,y))

  i = 1
  for (i in 1:length(id.r1)) {
    id = id.r1[i]
    id1 = edg %>% dplyr::filter( OTU_2 %in% c(id)) %>% .$OTU_1 %>% unique()
    id2 = edg %>% dplyr::filter( OTU_1 %in% c(id)) %>% .$OTU_2 %>% unique()
    id0 = c(id1,id2)

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
      dat1 = sz.node (x.0 = tem[1,1],y.0 = tem[1,2],r = r,n = n)
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

  dat0$class2 = "leaf"
  dat2$class2 = "brach"
  dat3 = dat0 %>% full_join(dat2)
  return(dat3)
}










square_points <- function(num_points, spacing) {
  # 计算正方形的边长
  side_length <- ceiling(sqrt(num_points))

  # 计算所有点的坐标
  y_coords <- rep(seq(0, (side_length - 1) * spacing, by = spacing), each = side_length)
  x_coords <- rep(seq(0, (side_length - 1) * spacing, by = spacing), times = side_length)[1:num_points]
  y_coords = sort(y_coords, decreasing = TRUE)[1:num_points]
  dat = data.frame(x =x_coords,y = y_coords )
  return(dat)
}


sz.node = function(
    x.0 = 0,
    y.0 = 0,
    r = 1,
    n = 46
){
  data = cbind(sin(2 * pi * ((0:(n - 1))/n))*r +x.0, cos(2 * pi * ((0:(n - 1))/n))*r +y.0)
  data =as.data.frame(data)
  colnames(data) = c("x","y")
  return(data)
}


