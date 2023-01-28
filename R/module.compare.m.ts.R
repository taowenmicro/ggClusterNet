

# res = module.compare.m.ts(
#   ps.st = ps.st,
#   N = 200,
#   degree = TRUE,
#   zipi = FALSE,
#   r.threshold= 0.8,
#   p.threshold=0.05,
#   method = "spearman",
#   padj = F,
#   n = 3,
#   g1 = "Group",# 分组1
#   g2 = "space",# 分组2
#   g3 = "time",# 分组3
#   zoom = 0.3,# 控制小圈大小
#   b.x = 1.3,
#   b.y = 1.3)
#
#
# p = res[[1]]
# #提取数据
# dat = res[[2]]
# #-作图数据提取
# dat = res[[3]]


module.compare.m.ts = function(
    ps.st = ps.st,
    N = 200,
    degree = TRUE,
    zipi = FALSE,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman",
    padj = F,
    n = 3,
    g1 = "Group",# 分组1
    g2 = "space",# 分组2
    g3 = "time",# 分组3
    ord.g1 =NULL,# 排序顺序
    ord.g2 = NULL, # 排序顺序
    ord.g3 = NULL,# 排序顺序
    zoom = 0.3,# 控制小圈大小
    b.x = 1,
    b.y = 1
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

      result = cor_Big_micro(ps = ps.g,
                             N = N,
                             r.threshold= r.threshold,
                             p.threshold= p.threshold,
                             method = method,
                             scale = FALSE)
      cor = result[[1]]
      tem = paste(dat.f[j,1],dat.f[j,2],group[n],sep = ".")
      cor.all[[tem]] = cor
      #--计算模块信息，部分OTU没有模块，注意去除
      result2 = model_maptree2(cor = cor,
                               method = "cluster_fast_greedy"
      )

      mod1 = result2[[2]]
      head(mod1)
      tem = paste(dat.f[j,1],dat.f[j,2],group[n],sep = ".")
      mod1 = mod1 %>% filter(!group == "mother_no") %>% select(ID,group)
      mod1$group = paste(tem,mod1$group,sep = "")
      mod1$Group = tem

      if (j ==1 & n == 1) {
        dat = mod1
      } else{
        dat = rbind(dat,mod1)
      }


    }
  }

  head(dat)

  node_table2  = dat
  head(node_table2)
  node_table2$Group %>% table()


  dat2 = model_compare(
    node_table2 = dat,
    n = n,
    padj = padj
  )


  # head(dat2)
  # head(node_table2)

  tem = node_table2 %>% distinct( group, .keep_all = TRUE)
  edge = data.frame(from = dat2$module1,to = dat2$module2,Value = 1)
  head(edge)
  id = c(tem$group) %>% unique()
  cor = matrix(0,nrow = length(id),ncol = length(id))
  colnames(cor) = id
  row.names(cor) = id

  netClu = data.frame(ID = tem$group,group = tem$Group)
  head(netClu)
  netClu$ID %>%unique()


  result2 = model_filled_circle(cor = cor,culxy =TRUE,
                                da = NULL,# 数据框，包含x,和y列
                                nodeGroup = netClu,
                                seed = 10,
                                mi.size = 0.2,
                                zoom = zoom)
  node = result2[[1]]
  head(node)
  head(edge)

  branch = result2[[2]]

  edge2 = edge %>% left_join(node,by = c("from" = "elements")) %>%
    dplyr::rename(x1 = X1,y1 = X2) %>%left_join(node,by = c("to" = "elements")) %>%
    dplyr::rename(x2 = X1,y2 = X2)
  head(edge2)

  ### 出图
  # library(ggrepel)

  pnet <- ggplot() +
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
                 data = edge2, size = 0.5,color = "#FF7F00") +
    geom_point(aes(X1, X2),pch = 21, data = node,fill ="#984EA3" ) +
    scale_colour_brewer(palette = "Set1") +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    # labs( title = paste(layout,"network",sep = "_"))+
    # ggrepel::geom_text_repel(aes(X1, X2,label=elements),size=1, data = node)+
    ggrepel::geom_text_repel(aes(x*b.x, y*b.y,label=group),size=4, data = branch)+
    theme(panel.background = element_blank()) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  pnet

  dat.all = list(nodule.all.g = node_table2,
                 compare.m.g = dat2
                 )
  plot.data = list(
    node = node,
    branch = branch,
    edge = edge2
  )
  return(list(pnet,dat.all,plot.data))
}
