

# library(phyloseq)
# library(ggClusterNet)
# library(tidyverse)
# library(ggnewscale)
# library(ggrepel)


# res = rank.network(
#     ps.st= ps.st,# phyloseq对象
#     g1 = "Group",# 分组1
#     g2 = "space",# 分组2
#     g3 = "time",# 分组3
#     ord.g1 = NULL, # 排序顺序
#     ord.g2 = NULL, # 排序顺序
#     ord.g3 = NULL, # 排序顺序
#     order = "space", # 出图每行代表的变量
#     jj = "Phylum",
#     fill = "Phylum",
#     method = "spearman",
#     clu_method = "cluster_fast_greedy",
#     select_layout = TRUE,
#     r.threshold=0.8,
#     p.threshold=0.01,
#     N= 500)
#
# p = res[[1]]
# ggsave("cs1.pdf",p,width = 5*9,height = 4*2)


rank.network = function(
    ps.st = ps.st,# phyloseq对象
    g1 = "Group",# 分组1
    g2 = "space",# 分组2
    g3 = "time",# 分组3
    ord.g1 = c("WT","KO","OE"), # 排序顺序
    ord.g2 = NULL, # 排序顺序
    ord.g3 = c("T1","T2","T3"), # 排序顺序
    order = "space", # 出图每行代表的变量
    jj = "Phylum",
    fill = "Phylum",
    method = "spearman",
    clu_method = "cluster_fast_greedy",
    select_layout = TRUE,
    r.threshold=0.8,
    p.threshold=0.01,
    N= 500

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
  cor.all2 = list()

  j = 1
  n = 2

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
      # head(map)
      # map$Group
      ps.g = ps.f  %>% subset_samples.wt(g1,group[n])

      result = cor_Big_micro(ps = ps.g,
                             N = N,
                             r.threshold= r.threshold,
                             p.threshold= p.threshold,
                             method = method,
                             scale = FALSE)
      cor = result[[1]]
      tem = paste(dat.f[j,1],dat.f[j,2],group[n],sep = ".")
      cor.all[[tem]] = cor


      otu_table = as.data.frame(t(vegan_otu(ps.g)))
      tax_table = as.data.frame(vegan_tax(ps.g))
      netClu = data.frame(ID = row.names(tax_table),group = tax_table[[jj]] )
      netClu$group = as.factor(netClu$group)
      result2 = model_maptree_group(
        cor = cor,
        nodeGroup =netClu,
        seed = 12)
      node = result2[[1]]

      #-----计算边
      edge = edgeBuild(cor = cor,node = node)
      netClu$group %>% unique()

      #-去除了分类等级内部相关
      #-内部OTU之间的相关都去除
      for (i in 1: length(levels(netClu$group))) {
        tem1 = levels(netClu$group)[i]
        tem2 = netClu$ID[netClu$group == tem1]

        if (i == 1) {
          tem3 = edge %>%
            filter(!(OTU_2 %in% tem2 & OTU_1 %in% tem2))
          dim(tem3)
        } else if(i > 1) {
          tem3 = tem3 %>%
            filter(!(OTU_2 %in% tem2 & OTU_1 %in% tem2))
          dim(tem3)

        }
      }

      tem4 = tem3 %>% left_join(netClu,by = c("OTU_2" = "ID")) %>%
        rename(OTU_2g = group) %>%
        left_join(netClu,by = c("OTU_1" = "ID")) %>%
        rename(OTU_1g = group)
      head(tem4)

      tem4$tem = paste(tem4$OTU_1g,tem4$OTU_2g,sep = "@")
      tem4$weight = abs(tem4$weight)
      tem5 = tem4 %>% group_by(tem) %>%
        summarise(sum(weight)) %>%
        rename( weight = `sum(weight)`) %>%
        as.data.frame()
      tem5$from = sapply(strsplit(tem5$tem, "[@]"), `[`, 1)
      tem5$to = sapply(strsplit(tem5$tem, "[@]"), `[`, 2)
      tem5 = tem5 %>% select(from,to,weight)

      # 这里将weight的进行标准化，这个也已经不是传统意义
      tem5$weight = tem5$weight/sum(tem5$weight)
      #--列表边矩阵
      mat = tidyfst::df_mat(tem5,from,to,weight)
      mat[is.na(mat)] = 0
      #-重新做好了 cor矩阵
      # cor = mat
      cor.all2[[tem]] = mat

      ps_net = ps.g %>% tax_glom_wt(jj)
      ps_net = ps_net %>%
        subset_taxa.wt("OTU",colnames(mat))


      otu_table = as.data.frame(t(vegan_otu(ps_net)))
      tax_table = as.data.frame(vegan_tax(ps_net))
      #-模拟一个分组
      netClu = data.frame(ID = row.names(tax_table),group = "A")
      netClu$group = as.factor(netClu$group)

      set.seed(12)
      # library(sna)

      # --随机模块化布局-模块多的话，这个布局浪费是非常多的，因为要将每一个模块分开，迭代很多次
      result2 = PolygonClusterG(cor = mat,nodeGroup = netClu )
      nod = result2[[1]]
      head(nod)
      # ---node节点注释
      nod2 = nodeadd(plotcord =nod,otu_table = otu_table,tax_table = tax_table)
      head(nod2)
      #-----计算边
      edg = edgeBuild(cor = mat,node = nod)
      # head(edg)
      # head(edge)
      nod2$group = tem
      edg$group = tem

      if (j ==1& n == 1) {
        node.r = nod2
        edge.r = edg
      } else{
        node.r = rbind(nod2,node.r)
        edge.r = rbind(edg,edge.r)
      }

    }
  }


  head(edge.r)
  head(node.r)

  #-统计边的节点数量 node link
  tem = edge.r$group %>% table() %>% as.data.frame()
  colnames(tem) = c("group","links")
  i = 1
  id = edge.r$group %>% unique()
  aa = c()
  for (i in 1:length(id)) {
    aa[i] = edge.r %>% filter(group == id[i]) %>%
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

  node.r$group = factor(node.r$group,levels = a)
  edge.r$group = factor(edge.r$group,levels = a)
  tem3 = tem3[match(a,tem3$group),]
  tem3$label = factor(tem3$label,levels = tem3$label)
  edge.r = edge.r %>% left_join(tem3,by = "group")

  edge.r$label = factor(edge.r$label,levels = as.character(tem3$label))

  head(node.r)

  node.r = node.r %>% left_join(tem3,by = "group")
  node.r$label = factor(node.r$label,levels = as.character(tem3$label))



  net.dat = list(
    cortab = cor.all2,
    node = node.r,
    edge = edge.r
  )

  head(node.r)

  p1 <- ggplot() +
    geom_segment(aes(x = X1,
                     y = Y1,
                     xend = X2,
                     yend = Y2,
                     color = weight,
                     size = weight),
                 data = edge.r,alpha = 1) +
    scale_color_gradientn(colours =c("grey99","grey70"))+
    ggnewscale::new_scale_color() +
    geom_point(aes(X1, X2,fill = !!sym(fill)),pch = 21, data = node.r,size = 4) +
    geom_text(aes(X1, X2,label = elements),pch = 21, data = node.r) +
    facet_wrap(.~ label,scales="free_y",ncol = row.num ) +
    # geom_text(aes(X1, X2,label = elements),pch = 21, data = nodes) +
    scale_colour_manual(values = c("#377EB8","#E41A1C")) +
    # scale_size(range = c(2, 5)) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme(panel.background = element_blank()) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  p1

  return(list(network.plot = p1,network.data = net.dat))
}
