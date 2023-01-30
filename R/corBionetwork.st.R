

# ps.st= ps.merge# phyloseq对象
# g1 = "Group"# 分组1
# g2 = NULL# 分组2
# g3 = NULL# 分组3
# ord.g1 = NULL # 排序顺序
# ord.g2 = NULL # 排序顺序
# ord.g3 = NULL # 排序顺序
# order = NULL # 出图每行代表的变量
#
# fill = "filed"
# size = "igraph.degree"
# method = "spearman"
# clu_method = "cluster_fast_greedy"
# select_layout = TRUE
# layout_net = "model_maptree2"
# r.threshold=0.8
# p.threshold=0.01
# maxnode = 5
# N= 500
# scale = TRUE
# env = NULL
# bio = TRUE
# minsize = 4
# maxsize = 14
# lab = NULL
# label = TRUE
# data(psITS)
# #--细菌和真菌ps对象中的map文件要一样
# ps.merge <- ggClusterNet::merge16S_ITS(ps16s = ps16s,
#                                        psITS = psITS,
#                                        N16s = 500,
#                                        NITS = 500
# )
#
# ps.merge
# map =  phyloseq::sample_data(ps.merge)
# ps.merge %>% vegan_tax()
#
#
# res = corBionetwork.st(
#     ps.st= ps.merge,# phyloseq对象
#     g1 = "Group",# 分组1
#     g2 = NULL,# 分组2
#     g3 = NULL,# 分组3
#     ord.g1 = NULL, # 排序顺序
#     ord.g2 = NULL, # 排序顺序
#     ord.g3 = NULL, # 排序顺序
#     order = NULL, # 出图每行代表的变量
#     fill = "filed",
#     size = "igraph.degree",method = "spearman",
#     clu_method = "cluster_fast_greedy",
#     select_layout = TRUE,layout_net = "model_maptree2",
#     r.threshold=0.8,
#     p.threshold=0.01,
#     maxnode = 5,
#     N= 500,scale = TRUE,env = NULL,
#     bio = TRUE,minsize = 4,maxsize = 14)
#
# res[[1]]
# res[[2]]

corBionetwork.st = function(
    ps.st= ps.merge,# phyloseq对象
    g1 = "Group",# 分组1
    g2 = NULL,# 分组2
    g3 = NULL,# 分组3
    ord.g1 = NULL, # 排序顺序
    ord.g2 = NULL, # 排序顺序
    ord.g3 = NULL, # 排序顺序
    order = NULL, # 出图每行代表的变量
    fill = "filed",
    size = "igraph.degree",
    method = "spearman",
    lab = NULL,
    label = TRUE,
    clu_method = "cluster_fast_greedy",
    select_layout = TRUE,
    layout_net = "model_maptree2",
    r.threshold=0.8,
    p.threshold=0.01,
    maxnode = 5,
    N= 500,
    scale = TRUE,
    env = NULL,
    bio = TRUE,
    minsize = 4,
    maxsize = 14
){

  if (scale) {
    ps.st  = ps.st %>% ggClusterNet::scale_micro()
  }


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

  if (!is.null(env)) {
    colnames(env)[1] = "ID"
    env_sub <-  env[match(mapsub$ID,env$ID),]
    head(env_sub)
  }


  j = 1
  n = 1

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
      ps.g = ps.f  %>% subset_samples.wt(g1,group[n]) %>% filter_OTU_ps(N)
      big =- TRUE
      if (bio) {
        if (!is.null(env)) {

          if (big) {

            result <- corBiostripeBig(data =  env_sub,
                                      group = envGroup,
                                      ps = ps.g,
                                      r.threshold = r.threshold,
                                      p.threshold = p.threshold,
                                      method = method)
          } else {
            result <- corBiostripe(data =  env_sub,
                                   group = envGroup,
                                   ps = ps.g,
                                   r.threshold = r.threshold,
                                   p.threshold = p.threshold,
                                   method = method)
          }


          #-- extract cor matrix
          occor.r = result[[1]]
          tax = as.data.frame((ggClusterNet::vegan_tax(ps.g)))

          if (length(tax$filed) != 0) {
            group2 <- data.frame(SampleID = row.names(tax),Group = tax$filed)
          } else {
            group2 <- data.frame(SampleID = row.names(tax),Group = "OTU")
          }

          colnames(envGroup) <-c("SampleID","Group")
          netClu = rbind(envGroup,group2)
          colnames(netClu) <- c("ID","group")
        } else {

          if (big) {
            result <- corBiostripeBig(ps = ps.g,
                                      r.threshold = r.threshold, p.threshold = p.threshold, method = method)
          } else {
            result <- corBiostripe(ps = ps.g,
                                   r.threshold = r.threshold, p.threshold = p.threshold, method = method)
          }

          #-- extract cor matrix
          occor.r = result[[1]]
          tax = as.data.frame((ggClusterNet::vegan_tax(ps.g)))

          if (length(tax$filed) != 0) {
            group2 <- data.frame(SampleID = row.names(tax),Group = tax$filed)
          } else {
            group2 <- data.frame(SampleID = row.names(tax),Group = "OTU")
          }

          netClu = group2
          colnames(netClu) <- c("ID","group")
        }

      }

      head(netClu)

      cor = result[[1]]
      tem = paste(dat.f[j,1],dat.f[j,2],group[n],sep = ".")
      cor.all[[tem]] = cor
      # gru.all[[tem]] = cor

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

      netClu$group2 = tem
      head(nod)
      head(edg)


      if (j ==1 & n == 1) {
        node = nod
        edge = edg
        netClu2 = netClu
      } else{
        node = rbind(node,nod)
        edge = rbind(edg,edge)
        netClu2 = rbind(netClu,netClu2)
      }

    }

  }

#   head(edge)
#   head(node)
#  edge$group %>% unique()
#  netClu2$group2%>% unique()
# netClu2$group = as.factor(netClu2$group)

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


  if (is.null(order)){
    order = ""

  }

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
  } else{
    a = NULL
    row.num = length(group)
  }

  if (!is.null(a)) {
    node$group = factor(node$group,levels = a)
    edge$group = factor(edge$group,levels = a)
    tem3 = tem3[match(a,tem3$group),]
  }else{
    node$group = factor(node$group)
    edge$group = factor(edge$group)
  }

  head(edge)
  head(node)

  tax = ps.st %>% vegan_tax() %>% as.data.frame() %>% rownames_to_column("elements")
  node = node %>% left_join(tax,by = "elements")


  tem3$label = factor(tem3$label,levels = tem3$label)
  edge = edge %>% left_join(tem3,by = "group")
  head(edge)
  edge$label = factor(edge$label,levels = as.character(tem3$label))

  head(node)

  node = node %>% left_join(tem3,by = "group")
  node$label = factor(node$label,levels = as.character(tem3$label))



  net.dat = list(
    cortab = cor.all,
    node = node,
    edge = edge
  )



    p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                  data = edge, size = 0.3,alpha = 0.5) +
      geom_point(aes(x = X1, y = X2,size = !!sym(size),fill = !!sym(fill)),pch = 21, data =  node) +
      scale_colour_brewer(palette = "Set1") +
      scale_size(range = c(minsize, maxsize)) +
      scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = NULL) +
      # labs( title = paste(layout,"network",sep = "_")) +
      facet_wrap(.~ label,scales="free_y",ncol = row.num ) +
      theme_void()
    p0


    tem2 = node %>% dplyr::filter(elements %in% lab[[1]])
    if (label == TRUE ) {

      if (!is.null(lab)) {
        p1 <- p0 + ggrepel::geom_text_repel(aes(X1, X2,label= elements),size=4, data = tem2)
      } else {
        p1 <- p0 + ggrepel::geom_text_repel(aes(X1, X2,label= elements),size=4, data = node)
      }

    }

  return(list(network.plot = p0,network.data = net.dat))

}

