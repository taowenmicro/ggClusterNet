# library(ggClusterNet)
# library(igraph)
# library(ggClusterNet)
# library(tidyverse)
# library(phyloseq)
# data(ps)
# res  = module_display.2(
#   pst = ps,
#   r.threshold= 0.6,
#   p.threshold=0.05,
#   select.mod = c("model_1","model_2","model_3","model_4"),#选择指定模块可视化
#   Top = 500,
#   num = 5, # 模块包含OTU数量少于5个的不展示,
#   leg.col = 9
# )
# # 全部模块输出展示
# p1 = res[[1]]
# p1
#
# p2 = res[[2]]
# p2
#
# p2 = res[[3]]
# p2
#
#
# ggsave("cs.pdf",p2,width = 5,height = 8)
# dat = res[[4]]
# head(dat)
# # 输出模块边
# res$netdata.selt.mod$edge
# # 输出模块节点
# res$netdata.selt.mod$node

module_display.2 = function(
    pst = ps,
    Top = 500,
    r.threshold= 0.8,
    p.threshold=0.05,
    select.mod = c("model_1","model_2","model_3"),
    num = 5,
    leg.col = 4,
    method.clu = "cluster_walktrap"
){

  pst = pst %>%
    filter_taxa(function(x) sum(x ) > 0, TRUE) %>%
    scale_micro("rela") %>%
    filter_OTU_ps(Top)

  result = cor_Big_micro(ps = pst,
                         N = 0,
                         r.threshold= r.threshold,
                         p.threshold= p.threshold,
                         method = "spearman")

  cor = result[[1]]
  head(cor)

  #-计算模块
  result = model_maptree2(cor = cor, method =  method.clu)
  node = result[[1]]
  netClu = result[[2]]

  head(netClu)
  tem = netClu$group %>% table() %>% as.data.frame() %>%
    arrange(desc(Freq))
  colnames(tem)[1] = "model"
  head(tem)

  lab.t = tem %>% filter(model != "mother_no") %>% filter(Freq > num) %>%
    .$model %>% as.character()


  branch = result[[3]] %>% filter(!str_detect(elements,"other"))
  branch$elements = paste("model_",branch$elements,sep = "")
  head(branch)

 branch = branch %>% filter(elements %in% lab.t)



  # ---node节点注释
  nodes = nodeadd(plotcord =node,
                  otu_table = pst %>%
                    vegan_otu() %>%
                    t() %>%
                    as.data.frame(),
                  tax_table = pst %>% vegan_tax() %>%
                    as.data.frame())
  # head(nodes2)

  nodes2 = nodes %>% inner_join(netClu,by = c("elements" = "ID"))
  nodes2$group = paste("Model_",nodes2$group,sep = "")

  #-----计算边
  edge = edgeBuild(cor = cor,node = node)
  dim(edge)
  ### 出图
  pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                  data = edge, size = 0.5) +
    geom_point(aes(X1, X2,fill = group,size = mean),pch = 21, data = nodes2) +
    geom_text(aes(x = x, y = y,label = elements), data = branch) +
    scale_colour_brewer(palette = "Set1") +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    # labs( title = paste(layout,"network",sep = "_"))+
    # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
    # discard default grid + titles in ggplot2
    theme(panel.background = element_blank()) +
    # theme(legend.position = "none") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(legend.position = 'none')


  #--挑选模块展示
  mod1 = netClu
  # head(mod1)
  tem = mod1$group %>% table() %>%
    as.data.frame() %>%
    dplyr::arrange(desc(Freq))
  colnames(tem) = c("Model","OTU.num")
  head(tem)


  if (length(select.mod) == 1 & is.numeric(select.mod)) {
    select.mod.name = tem$Model[1:select.mod]
    mod1 = mod1 %>% filter(!group == "mother_no",
                           group %in%c(select.mod.name)

    ) %>% select(ID,group,degree)

  } else if (is.character(select.mod)) {
    select.mod.name = select.mod
    mod1 = mod1 %>% filter(!group == "mother_no",
                           group %in%c(select.mod.name)

    ) %>% select(ID,group,degree)

  }

  # head(mod1)
  # head(node)
  node = result[[1]] %>% filter(elements %in% mod1$ID)
  # ---node节点注释
  nodes = nodeadd(plotcord =node,
                  otu_table = pst %>%
                    vegan_otu() %>%
                    t() %>%
                    as.data.frame(),
                  tax_table = pst %>% vegan_tax() %>%
                    as.data.frame())
  head(nodes)

  nodes2 = nodes %>% inner_join(mod1,by = c("elements" = "ID"))

  #-----计算边
  edge = edgeBuild(cor = cor[mod1$ID,mod1$ID],node = node)

  ### 出图
  p2 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 0.5) +
    geom_point(aes(X1, X2,fill = group,size = mean),pch = 21, data = nodes2) +
    scale_colour_brewer(palette = "Set1") +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    theme_void()



  net.s = netClu[netClu$group %in% lab.t,]

  #-计算模块
  result = model_maptree_group(cor = cor, nodeGroup = net.s,
                               seed = 12)
  node = result[[1]]
  tem = net.s$group %>% table() %>% as.data.frame() %>%
    arrange(desc(Freq))
  colnames(tem)[1] = "model"


  branch = result[[2]] %>% filter(!str_detect(elements,"other"))
  branch$elements = paste("model_",branch$elements,sep = "")



  # ---node节点注释
  nodes = nodeadd(plotcord =node,
                  otu_table = pst %>%
                    vegan_otu() %>%
                    t() %>%
                    as.data.frame(),
                  tax_table = pst %>% vegan_tax() %>%
                    as.data.frame())
  # head(nodes2)

  nodes2 = nodes %>% inner_join(net.s,by = c("elements" = "ID"))
  # nodes2$group = paste("Model_",nodes2$group,sep = "")

  #-----计算边
  # cor[node$elements,node$elements]
  edge = edgeBuild(cor = cor[node$elements,node$elements],node = node)
  head(edge)

  head(net.s)
  n.1 = net.s$group %>% unique() %>% length()
  tab = data.frame(ID = net.s$group %>% unique(),color = RColorBrewer::brewer.pal(9,"Set1")[1:n.1])
  head(tab)
  tab$color

  tem = edge %>% left_join(net.s,by = c("OTU_1" = "ID")) %>%
    rename(group1 = group) %>%
    select(-degree) %>%
    left_join(net.s,by = c("OTU_2" = "ID")) %>%
    rename(group2 = group)
  head(tem)

  edge2 = tem %>% mutate(color1 = ifelse(group1 == group2,as.character(group1),"acorss")) %>%
    left_join(tab,by = c("color1" = "ID"))
  head(edge2)
  edge2$color[is.na(edge2$color)] = "grey80"
  # edge2$color = factor(edge2$color,levels = c(tab$color,"grey80"))
  # nodes2$group = factor(nodes2$group,levels = )


  edgeb = edge2 %>%
    filter(color1 == "acorss")
  head(edgeb)

  edge3 = edge2 %>%
    filter(color1 != "acorss")

  nodes2$group = as.factor(nodes2$group)
  edge3$color1 = factor(edge3$color1,levels = levels(nodes2$group))

  ### 出图
  p3 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2),color = edgeb$color,
                                data = edgeb, size = 0.5) +
    ggnewscale:: new_scale_color() +
    geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color1),
                   data = edge3, size = 0.5) +
    geom_point(aes(X1, X2,fill = group,size = mean),pch = 21, data = nodes2) +
    ggrepel::geom_text_repel(aes(x = x, y = y,label = elements), data = branch) +
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    theme_void() +
    theme(legend.position = "top") +
    guides(color=guide_legend(nrow=leg.col, byrow=TRUE),
           fill=guide_legend(nrow=leg.col, byrow=TRUE),
             size = guide_legend(nrow=3, byrow=TRUE),
           )

  p3

  netdata.selt.mod = list(edge = edge2,
                          node = nodes2
  )


  return(list(plot1 = pnet,plot2 = p2,plot3 = p3,
              cormatrix = cor,
              mod.groups = netClu,
              netdata.selt.mod = netdata.selt.mod))
}
