
# res = Facet.network (
#     ps.st= ps.st,# phyloseq对象
#     g1 = "Group",# 分组1
#     g2 = "space",# 分组2
#     g3 = "time",# 分组3
#     ord.g1 = c("WT","KO","OE"),# 排序顺序
#     ord.g2 = c("B","R") ,# 排序顺序
#     ord.g3 = c("T1","T2","T3") ,# 排序顺序
#     order = "space", # 出图每行代表的变量
#     fill = "Phylum",
#     size = "igraph.degree",
#     method = "spearman",
#     clu_method = "cluster_fast_greedy",
#     select_layout = TRUE,
#     layout_net = "model_maptree2",
#     r.threshold=0.8,
#     p.threshold=0.01,
#     maxnode = 5
# )
# p = res[[1]]
# nettab = res[[2]]

Facet.network = function(
    ps.st= ps.st,# phyloseq对象
    N = 200,
    g1 = "Group",# 分组1
    g2 = NULL,# 分组2
    g3 = NULL,# 分组3
    ord.g1 = NULL,# 排序顺序
    ord.g2 = NULL ,# 排序顺序
    ord.g3 = NULL ,# 排序顺序
    order = "space", # 出图每行代表的变量
    fill = "Phylum",
    size = "igraph.degree",
    method = "spearman",
    clu_method = "cluster_fast_greedy",
    select_layout = TRUE,
    layout_net = "model_maptree2",
    r.threshold=0.8,
    p.threshold=0.01,
    maxnode = 5,
    R = 100,
    ncpus = 1
  ){

  ps.all = ps.st
  map = sample_data(ps.all)

  # g2 = NULL
  if (is.null(g2)) {
    sp = ""
  } else if (is.null(ord.g2)){
    sp = map[,g2] %>% as.matrix() %>% as.vector() %>% unique()
  } else{
    sp = ord.g2
    print(sp)
  }
  # g3 = NULL
  if (is.null(g3)) {
    ti = ""
  } else if (is.null(ord.g3)){
    ti = map[,g3] %>% as.matrix() %>% as.vector() %>% unique()
  } else{
    ti = ord.g3
    print(ti)
  }

  if (is.null(ord.g1)) {
    group = map[,g1] %>% as.matrix() %>% as.vector() %>% unique()

  } else{
    group = ord.g1
    print(group)
  }



  #-构造两两组合全部情况
  for (i in 1:length(sp)) {
    dat = data.frame(g2 = sp[i],g3 = ti)

    if (i ==1) {
      dat.f = dat
    } else{
      dat.f = rbind(dat.f,dat)
    }

  }

  # head(dat.f)

  cor.all = list()
  # 存储不同分组的相关矩阵拮即可
  for (j in 1:nrow(dat.f)) {
    if (dat.f[j,1] == "") {
      ps.t = ps.all
    } else{
      ps.t = ps.all %>% subset_samples.wt(g2,dat.f[j,1])
    }

    if (dat.f[j,2] == "") {
      ps.f = ps.t
    } else{
      ps.f = ps.t  %>% subset_samples.wt(g3,dat.f[j,2])
    }


    for (n in 1:length(group)) {
      map = sample_data(ps.f)
      head(map)
      map$Group
      ps.g = ps.f  %>% subset_samples.wt(g1,group[n])

      if (method != "sparcc") {
        result = cor_Big_micro(ps = ps.g,
                               N = N,
                               r.threshold= r.threshold,
                               p.threshold= p.threshold,
                               method = method,
                               scale = FALSE)
        cor = result[[1]]
        print("cor matrix culculating over")
      } else if (method %in% c("sparcc")){
        result = corMicro (ps = psi,N = 0,r.threshold= r.threshold,p.threshold=p.threshold,
                           method = method,R = R,ncpus = ncpus)
        a = 1

      print("cor matrix culculating over")
      cor = result[[1]]

      }



      tem = paste(dat.f[j,1],dat.f[j,2],group[n],sep = ".")
      cor.all[[tem]] = cor

      #--构造边和节点文件，全部放到一起
      res = node.edge(
        cor = cor,
        select_layout = T,
        clu_method=clu_method,
        layout_net = layout_net
      )

      nod = res[[1]]
      nod$group = tem
      edg = res[[2]]
      edg$group = tem

      if (j ==1 & n == 1) {
        node = nod
        edge = edg
      } else{
        node = rbind(node,nod)
        edge = rbind(edg,edge)
      }

    }

  }



  tax = ps.all %>% vegan_tax() %>% as.data.frame() %>%
    rownames_to_column("ID")
  node$ID = node$elements
  node.1 = node  %>% left_join(tax,by = "ID")
  head(node.1)

  node.1$Group = sapply(strsplit(node.1$group, "[.]"), `[`, 3)
  node.1$time = sapply(strsplit(node.1$group, "[.]"), `[`, 2)
  node.1$space = sapply(strsplit(node.1$group, "[.]"), `[`, 1)

  edge$Group = sapply(strsplit(edge$group, "[.]"), `[`, 3)
  edge$time = sapply(strsplit(edge$group, "[.]"), `[`, 2)
  edge$space = sapply(strsplit(edge$group, "[.]"), `[`, 1)

  head(edge)
  #-统计边的节点数量 node link
  tem = edge$group %>% table() %>% as.data.frame()
  colnames(tem) = c("group","links")
  i = 1
  id = edge$group %>% unique()
  aa = c()
  for (i in 1:length(id)) {
    aa[i] = edge %>% filter(group == id[i]) %>%
      select("OTU_2", "OTU_1") %>% as.matrix() %>%
      as.vector() %>% unique() %>% length()
  }
  tem2 = data.frame(group = id,nodes = aa)

  tem3 = tem %>% full_join(tem2,by = "group")
  tem3$label= paste(tem3$group,": (nodes: ",
                    tem3$nodes,"; links: ",tem3$links,")",sep = "")


  #-行--空间
  if (order == "space"|order == "g2") {
    row.id = g3
    row.num = length(group) * length(ti)
    a = c()
    for (i in 1:length(sp)) {
      for (j in 1:length(ti)) {
        tem = paste(sp[i],ti[j],group,sep = ".")
        a = c(a,tem)
      }
    }

  } else if (order == "time"|order == "g3") {
    row.id = g2
    row.num = length(group) * length(sp)
    a = c()
    for (j in 1:length(ti)) {
      for (i in 1:length(sp)) {
        tem = paste(sp[i],ti[j],group,sep = ".")
        a = c(a,tem)
      }
    }
  }

  node.1$group = factor(node.1$group,levels = a)
  edge$group = factor(edge$group,levels = a)
  tem3 = tem3[match(a,tem3$group),]
  tem3$label = factor(tem3$label,levels = tem3$label)
  edge = edge %>% left_join(tem3,by = "group")
  head(edge)
  edge$label = factor(edge$label,levels = as.character(tem3$label))

  head(node.1)

  node.1 = node.1 %>% left_join(tem3,by = "group")
  node.1$label = factor(node.1$label,levels = as.character(tem3$label))



  net.dat = list(
    cortab = cor.all,
    node = node.1,
    edge = edge
  )

  # edge$group %>% unique()
  p <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,
                                   color = cor),
                               data = edge, size = 0.03,alpha = 0.5) +
    geom_point(aes(X1, X2,
                   fill = !!sym(fill),
                   size = !!sym(size) ),
               pch = 21, data = node.1,color = "gray40") +
    facet_wrap(.~ label,scales="free_y",ncol = row.num ) +
    # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
    # geom_text(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
    scale_colour_manual(values = c("#6D98B5","#D48852")) +
    scale_size(range = c(0.8, maxnode)) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5)
    ) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()
          ) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

  return(list(network.plot = p,network.data = net.dat))
}



#--这种方式似乎只能全部分组组合都要出图到一张上
# 因此需要人工可以指定排序方向和行列数量



node.edge = function(
  cor = cor,
  select_layout = T,
                     clu_method=clu_method,
                     layout_net = layout_net
){

  if (select_layout) {
    node = NULL
    node = culculate_node_axis(
      cor.matrix = cor,
      layout = layout_net,
      seed = 1,
      group = NULL,
      model = FALSE,
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





