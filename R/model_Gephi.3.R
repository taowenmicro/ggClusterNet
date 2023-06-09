# rm(list=ls())
# library(ggClusterNet)
# library(tidyverse)
# library(phyloseq)
# #----------计算相关#----
# result = cor_Big_micro(ps = ps,
#                        N = 2000,
#                        # method.scale = "TMM",
#                        r.threshold=0.5,
#                        p.threshold=0.05,
#                        method = "spearman"
# )
#
# # result = corMicro (ps = ps,
# #                    N = 500,
# #                    r.threshold= 0.8,
# #                    p.threshold=0.05,
# #                    method = "sparcc",
# #                    R = 10,
# #                    ncpus = 5)
#
#
# #--提取相关矩阵
# cor = result[[1]]
# library(igraph)
#
# plots = list()
# id = c(0,0.2,0.4,0.6,0.8,1)
# for (j in 1:length(id)) {
#   i = id[j]
#   # tab = model_Gephi.3(
#   #   cor = cor,
#   #   t0 = 0.98, # 取值范围为0到1，表示位于两个已知点之间的中间位置
#   #   t2 = 1.2, # 每个聚类中心的缩放
#   #   t3 = 8) # 聚类取哪个点，数值越大，取点越靠近00点
#
#
#   tem3 = tab[[1]]
#   node.p = tab[[2]]
#   head(tem3)
#
#
#   aa = node.p$group %>% table() %>% as.data.frame() %>%
#     rename(  "group" = ".")
#
#   tid = aa$group[aa$Freq < 60] %>% as.character()
#   tem3$group = as.character(tem3$group)
#   tem3$group[(tem3$group) %in% tid] = "miniModule"
#   test_color5 <- RColorBrewer::brewer.pal(n = 12,name = "Set3")
#   head(tem3)
#   tem3$ID = NULL
#   colnames(tem3)[5] = "Group"
#   head(node.p)
#
#   nodef = node.p %>% left_join(tem3,by = c("ID" = "OTU"))
#   head(nodef)
#
#   # p1 <- ggplot() +
#   #   # geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2),
#   #   #                             data = edge, size = 1) +
#   #   geom_point(aes(X1, X2,fill = group), data = tem3,pch = 21,
#   #              position=position_jitter(width=0.15,height=0.15)
#   #
#   #              ) +
#   #   theme_void()
#   # p1
#
#   # "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF" "#999999"
#   p2 <- ggplot() +
#     # geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
#     #                               data = edge, size = 0.5,alpha = 0.01) +
#     geom_point(aes(X1, X2,fill = Group,size = igraph.degree),pch = 21, data = nodef,color= "grey80",
#                position=position_jitter(width=i,height=i)
#     ) +
#     scale_colour_brewer(palette = "Set1") +
#     scale_fill_hue() +
#     scale_fill_manual(values  =
#                         c(
#                           "#4DAF4A","#A65628","#FDAE61" ,"#FEE08B" ,"#3288BD","grey80","grey90")
#     )+
#     scale_x_continuous(breaks = NULL) +
#     scale_y_continuous(breaks = NULL) +
#     scale_size(range = c(3, 20)) +
#     # labs( title = paste(layout,"network",sep = "_"))+
#     # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
#     # discard default grid + titles in ggplot2
#     theme(panel.background = element_blank()) +
#     # theme(legend.position = "none") +
#     theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
#     theme(legend.background = element_rect(colour = NA)) +
#     theme(panel.background = element_rect(fill = "white",  colour = NA)) +
#     theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
#   p2
#   plots[[j]] = p2
# }
#
# # ggsave("cs24.pdf",p2,width = 20,height = 19)
#
#
#
# p1  = ggpubr::ggarrange(plotlist = plots,
#                         common.legend = FALSE, legend="right",ncol = 3,nrow = 2)
# p1
# ggsave("cs2.pdf",p1,width = 20*3,height = 19*2,limitsize = FALSE)



model_Gephi.3 = function(
    cor = cor,
    t0 = 1, # 取值范围为0到1，表示位于两个已知点之间的中间位置
    t2 = 1.2, # 每个聚类中心的缩放
    t3 = 4 # 聚类取哪个点，数值越大，取点越靠近00点
){

  # 给出坐标数据
  num.node <- dim(cor)[1]
  for (N in 1: num.node) {
    A = 1 + (7*(N + 1)*N )/2 - N
    if (A >= num.node) {

      break
    }
    n = N - 1
    # print(n)
  }


  # n = (sqrt((num.node-1)/3) - 1) %>% floor()
  wai.mode = num.node - (1 + (7*(n + 1)*n )/2 - n)
  dat = data.frame(x = 0,y = 0)
  for (i in 1:n) {
    t <- seq(0, 2*pi, length.out = 7*i)
    t = t[-1]
    x <- sin(t)*i
    y <- cos(t)*i
    add = data.frame(x = x,y = y)
    dat = rbind(dat,add)

    if (i== n) {
      i = i + 1
      t <- seq(0, 2*pi, length.out = (wai.mode + 1))
      t = t[-1]
      x <- sin(t)*i
      y <- cos(t)*i
      add = data.frame(x = x,y = y)
      dat = rbind(dat,add)
    }

  }
  row.names(dat) = row.names(cor)
  dat$elements = row.names(cor)
  colnames(dat)[1:2] = c("X1","X2")
  head(dat)

  #--给出填充大小的标度加入算法，使得大小相似的点位于距离较为近的地方

  result2 = model_Gephi.2(cor = cor)
  node = result2[[3]]
  head(node)
  igraph = make_igraph(cor)

  node.p = node_properties(igraph) %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    full_join(node,by = "ID")
  head(node.p)
  node.p$igraph.degree[is.na(node.p$igraph.degree)] = 0.1


  clutab = dat[,1:2]
  ztab = data.frame(ID = node.p$ID,size = node.p$igraph.degree,group = node.p$group)
  head(ztab)

  # tab3 = node.p$group %>% table() %>% as.data.frame() %>%
  #   rename(  "group" = ".") %>%
  #   arrange(desc(Freq))
  node.p$tem = 1
  tab3 = node.p %>% group_by(group) %>%
    summarise(Freq = sum(tem),numedge = sum(igraph.degree)) %>%
    arrange(desc(numedge))
  head(tab3)

  # tab3$radio = tab3$Freq/tab3$numedge
  #   tab3 = tab3 %>%  arrange(desc(radio))


  # #--第一个布局#---------
  # t <- 1  # 取值范围为0到1，表示位于两个已知点之间的中间位置
  # t2 = 1.2 # 每个聚类中心的缩放
  # t3 = 4
  # i= 1
  for (i in 1:length(tab3$group)) {
    # 修改一下标签，防止和微生物名字冲突
    if (i == 1) {
      row.names(clutab) = paste0("A",1:length(row.names(clutab)))
      xytab = clutab[,1:2]
      row.names(xytab) = paste0("A",1:length(row.names(xytab)))
    }
    # # 使用k-means算法进行聚类
    # if (dim(clutab)[1]< 4 ) {
    #
    # }
    # 计算模块位置
    if (dim(clutab)[1] > (length(tab3$group) +2 -i)) {
      kmeans_result <- kmeans(clutab, centers = (length(tab3$group) +2-i), nstart = 20, iter.max = 100)
      cluster_centers <- kmeans_result$centers
      aa = apply(kmeans_result$centers, 1, function(x) sqrt(sum((x - cluster_centers)^2))) %>% as.data.frame()
      colnames(aa) = "dis"
      aa = aa %>% rownames_to_column("id") %>%
        arrange(dis)
      center = kmeans_result$centers[aa$id,]

      tem = c(0,0)
      tem1 = dist(rbind(tem,center)) %>% as.matrix() %>% tidyfst::mat_df() %>%filter(row == "tem") %>%
        arrange(value)

      center = center[-as.numeric(tem1[2,2]),]
    } else if(dim(clutab)[1] == (length(tab3$group) +2-i)) {

      kmeans_result <- kmeans(clutab, centers = (length(tab3$group) +1-i), nstart = 20, iter.max = 100)
      cluster_centers <- kmeans_result$centers
      aa = apply(kmeans_result$centers, 1, function(x) sqrt(sum((x - cluster_centers)^2))) %>% as.data.frame()
      colnames(aa) = "dis"
      aa = aa %>% rownames_to_column("id") %>%
        arrange(dis)
      center = kmeans_result$centers[aa$id,]

      tem = c(0,0)
      tem1 = dist(rbind(tem,center)) %>% as.matrix() %>% tidyfst::mat_df() %>%filter(row == "tem") %>%
        arrange(value)

      center = center
    }

    if (i == length(tab3$group)) {
      cluster_centers <- center
    } else if(i %in% c(1:3)){
      aa = apply(center, 1, function(x) sqrt(sum((x - c(0,0))^2))) %>% as.data.frame() %>%
        rename(  "dis" = ".")
      aa$id = row.names(center)

      if (t3 == 0) {
        a2 = 1
      } else {
        a2 = nrow(center)%/%t3
      }

      if (a2 == 0) {
        a2 = 1
      }


      cluster_centers <- center[aa %>% arrange((dis)) %>% .$id %>% .[a2],]


    } else{
      aa = apply(center, 1, function(x) sqrt(sum((x - c(0,0))^2))) %>% as.data.frame() %>%
        rename(  "dis" = ".")
      aa$id = row.names(center)
      cluster_centers <- center[aa %>% arrange((dis)) %>% .$id %>% .[1],]
    }
    # 改变中心点到整个网络中心的距离
    x2 = cluster_centers[1]
    y2 = cluster_centers[2]
    x1 = 0
    y1 = 0
    # 定义所需点的位置参数
    cluster_centers[1] <- (1 - t2) * x1 + t2 * x2
    cluster_centers[2] <- (1 - t2) * y1 + t2 * y2

    # 指定要聚为一类的数量
    num_points_in_cluster <- as.data.frame(tab3)[i,2]
    distances <- apply(clutab, 1, function(x) sqrt(sum((x - cluster_centers)^2))) %>% as.data.frame()
    colnames(distances)[1] = "dis"

    head(distances)
    id = distances %>%
      dplyr::arrange(dis) %>%
      dplyr::slice(1:num_points_in_cluster)%>%
      row.names()

    tem = xytab %>% rownames_to_column("ID")
    tem2 = distances %>%
      rownames_to_column("ID" ) %>%
      dplyr::arrange(dis) %>%
      dplyr::slice(1:num_points_in_cluster)  %>%
      left_join(tem,by = "ID") %>%
      dplyr::arrange((dis))
    head(tem2)
    tem2$OTU = ztab %>% filter(group == as.data.frame(tab3)[i,1] ) %>% dplyr::arrange(desc(size)) %>%
      .$ID
    tem2$group = as.data.frame(tab3)[i,1]

    head(tem2)
    # 这里缩放每个聚类模块相对大小
    x1 = cluster_centers[1]
    y1 = cluster_centers[2]
    # j = 1
    # t <- 0.95  # 取值范围为0到1，表示位于两个已知点之间的中间位置
    for (j in 1:nrow(tem2)) {
      x2 = tem2[j,3]
      y2 = tem2[j,4]
      # 定义所需点的位置参数
      tem2[j,3] <- (1 - t0) * x1 + t0 * x2
      tem2[j,4] <- (1 - t0) * y1 + t0 * y2
    }



    clutab= clutab %>% rownames_to_column("ID") %>%
      filter(!ID %in% id) %>% column_to_rownames("ID")

    if (i ==1) {
      tem3 = tem2
    } else{
      tem3 = rbind(tem3,tem2)
    }
  }
  return(list(tem3,node.p))
}

