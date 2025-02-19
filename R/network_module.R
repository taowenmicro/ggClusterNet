
# res = net.property.module.env (
#     pst = pst,
#     Top = 500,
#     r.threshold= 0.8,
#     p.threshold=0.05,
#     env = env,
#     select.mod = c("model_1","model_2","model_3"),
#     select.env = "pH")
#
# res[[1]]
# res[[2]]
# res[[3]]

net.property.module.env = function(
    pst = ps,
    corg = NULL,
    method = "pearson",
    Top = 500,
    r.threshold= 0.8,
    p.threshold=0.05,
    env = NULL,
    select.mod = c("model_1","model_2","model_3"),
    select.env = NULL
){

  # result = cor_Big_micro(ps = pst,
  #                        N = 0,
  #                        r.threshold= r.threshold,
  #                        p.threshold= p.threshold,
  #                        method = "spearman")
  #
  # cor = result[[1]]
  # head(cor)

  if (is.null(corg)) {
    # pst =  ps %>%
    #   scale_micro() %>%
    #   subset_samples.wt("Group", c(id[i])) %>%
    #   filter_OTU_ps(Top)

    result = cor_Big_micro(ps = pst,
                           N = Top,
                           r.threshold= r.threshold,
                           p.threshold= p.threshold,
                           method = method)

    cor = result[[1]]
    # head(cor)
  } else if (!is.null(corg)){
    cor = corg
  }



  result2 = model_maptree2(cor = cor, method = "cluster_fast_greedy")

  # select.mod = select.mod
  mod1 = result2[[2]]
  # head(mod1)
  tem = mod1$group %>% table() %>%
    as.data.frame() %>%
    dplyr::arrange(desc(Freq))
  colnames(tem) = c("Model","OTU.num")

  # head(tem)
  if (length(select.mod) == 1 & is.numeric(select.mod)) {
    select.mod.name = tem$Model[1:select.mod]
    mod1 = mod1 %>% filter(!group == "mother_no",
                           group %in%c(select.mod.name)

    ) %>% select(ID,group,degree)

  } else if (is.character(select.mod)& select.mod[1] != "no") {
    select.mod.name = select.mod
    mod1 = mod1 %>% filter(!group == "mother_no",
                           group %in%c(select.mod.name)

    ) %>% select(ID,group,degree)

  }else if (select.mod == "no") {
    select.mod.name = select.mod
    mod1 = mod1 %>% filter(!group == "mother_no")

  }

  id.s = mod1$group %>% unique()
  for (i in 1:length(id.s)) {
    id.t =  mod1 %>%
      dplyr::filter(group %in% id.s[i]) %>%
      .$ID
    ps.t = ps %>%
      scale_micro() %>%
      subset_taxa.wt("OTU", id.t )

    otu = ps.t %>%
      vegan_otu() %>%
      t()

    colSD = function(x){
      apply(x,2, sd)
    }

    dat = (otu - colMeans(otu))/colSD(otu)
    head(dat)
    otu_table(ps.t) = otu_table(as.matrix(dat),taxa_are_rows = T)
    #--计算总丰度
    otu = ps.t %>%  vegan_otu() %>% t()

    colSums(otu)

    dat = data.frame(id = names(colSums(otu)),abundance.zscore = colSums(otu))
    colnames(dat)[2] = id.s[i]

    if (i ==1) {
      tem = dat
    } else{
      dat$id = NULL
      tem = cbind(tem,dat)
    }
  }

  head(tem)
  map =sample_data(ps.t)
  map$id = row.names(map)
  map = map[,c("id","Group")]
  data = map %>%
    as.tibble() %>%
    inner_join(tem,by = "id") %>%
    dplyr::rename(group = Group)

  if (!is.null(env)) {
    colnames(env)[1] = "id"
    subenv = env %>% dplyr::select(id,everything()) %>% dplyr::select(id,select.env )
    # head(data)
    tab = data %>% left_join(subenv,by = "id")
    modenv = tab
    # head(tab)
    # library(reshape2)
    mtcars2 = reshape2::melt(tab, id.vars=c(select.env,"group","id"))
    mtcars2$variable
    head(mtcars2)
    lab = mean(mtcars2[,select.env])
    p1_1 = ggplot2::ggplot(mtcars2,aes(x= value,!!sym(select.env), colour=variable)) +
      ggplot2::geom_point() +
      ggpubr::stat_cor(label.y=lab*1.1)+
      ggpubr::stat_regline_equation(label.y=lab*1.1,vjust = 2) +
      facet_wrap(~variable, scales="free_x") +
      geom_smooth(aes(value,!!sym(select.env), colour=variable), method=lm, se=T)+
      theme_classic()
  } else{
    modenv =  NULL
    p1_1 = NULL
  }

  # p1_1

  dat.f = netproperties.sample(pst = pst,cor = cor)
  # head(dat.f)
  dat.f$id = row.names(dat.f)
  dat.f = dat.f %>% dplyr:: select(id,everything())
  tab = dat.f %>% left_join(subenv,by = "id")
  # head(tab)

 if (!is.null(env)) {
   mtcars2 = reshape2::melt(tab, id.vars=c(select.env,"id"))
   lab = mean(mtcars2[,select.env])
   # head(mtcars2)
   preoptab = tab

   p0_1 = ggplot2::ggplot(mtcars2,aes(x= value,!!sym(select.env), colour=variable)) +
     ggplot2::geom_point() +
     ggpubr::stat_cor(label.y=lab*1.1)+
     ggpubr::stat_regline_equation(label.y=lab*1.1,vjust = 2) +
     facet_wrap(~variable, scales="free_x") +
     geom_smooth(aes(value,!!sym(select.env), colour=variable), method=lm, se=T)+
     theme_classic()
 } else{
   preoptab = NULL
   p0_1 = NULL
 }

  # p0_1
  plotdat = list(
    model.env = modenv,
    preopertites.env = preoptab

  )
  return(list(p1_1,p0_1,dat.f,plotdat))
}



# dat.f = netproperties.sample(pst = pst,cor = cor)

netproperties.sample = function(
    pst = pst,
    cor = cor){
  igraph = make_igraph(cor)
  dat = igraph::V(igraph)
  names(dat) %>% length()
  otu = pst %>% vegan_otu() %>% t()
  otu = otu[row.names(otu) %in% names(dat),]


  otu[otu > 1] = 1
  dim(otu)
  A = list()
  dat.f = NULL



  for (i in 1:length(colnames(otu))) {
    tem = otu[,colnames(otu)[i]][otu[,colnames(otu)[i]] > 0 ] %>% names()
    A[[colnames(otu)[i]]] = tem
    #-计算性质
    tem.2 = A[[colnames(otu)[i]]]
    tem.g = igraph::induced_subgraph(igraph,tem.2)
    dat = net_properties.2(tem.g,n.hub = FALSE)
    head(dat,n = 16)

    dat[16,1] = 0
    dat = as.data.frame(dat)
    dat$value = as.numeric(dat$value)
    colnames(dat) = colnames(otu)[i]
    if (i == 1) {
      dat.f = dat
    } else {
      dat.f = cbind(dat.f,dat)
    }
  }
  head(dat.f)

  dat.f = dat.f %>%
    t() %>%
    as.data.frame()
  return(dat.f)
}










# res = module_abundance(
#     pst = pst,
#     mod1 = mod1
# )
# res[[1]]
# res[[2]]


module_abundance = function(
    pst = pst,
    mod1 = mod1
){
  id.s = mod1$group %>% unique()

  for (i in 1:length(id.s)) {
    id.t =  mod1 %>%
      dplyr::filter(group %in% id.s[i]) %>%
      .$ID
    ps.t = ps %>%
      scale_micro() %>%
      subset_taxa.wt("OTU",id.t )

    otu = ps.t %>%
      vegan_otu() %>%
      t()



    colSD = function(x){
      apply(x,2, sd)
    }

    dat = (otu - colMeans(otu))/colSD(otu)
    head(dat)
    otu_table(ps.t) = otu_table(as.matrix(dat),taxa_are_rows = T)

    #--计算总丰度

    otu = ps.t %>%  vegan_otu() %>% t()

    colSums(otu)

    dat = data.frame(id = names(colSums(otu)),abundance.zscore = colSums(otu))
    colnames(dat)[2] = id.s[i]

    if (i ==1) {
      tem = dat
    } else{
      dat$id = NULL
      tem = cbind(tem,dat)
    }
  }

  head(tem)
  map =sample_data(ps.t)
  map$id = row.names(map)
  map = map[,c("id","Group")]
  data = map %>%
    as.tibble() %>%
    inner_join(tem,by = "id") %>%
    dplyr::rename(group = Group)

  num = 3+ length(id.s) -1
  result = EasyStat::MuiaovMcomper2(data = data,num = c(3:num))


  result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = c(3:num),
                                            result = result,
                                            sig_show ="abc",ncol = length(id.s) )
  p1_1 = result1[[1]] +
    ggplot2::guides(fill = guide_legend(title = NULL)) +
    theme_classic()
  # p1_1
  p1_2 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
    geom_violin(alpha=1, aes(fill=group)) +
    geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
    labs(x="", y="Module Abundance (z-score)")+
    facet_wrap(.~name,scales="free_y",ncol  = length(id.s)) +
    # theme_classic()+
    geom_text(aes(x=group , y=y ,label=stat)) +

    guides(color=guide_legend(title = NULL),
           shape=guide_legend(title = NULL),
           fill = guide_legend(title = NULL)
    ) +
    theme_classic()
  # p1_2
  return(list(p1_1,p1_2,result1[[2]]))
}




# res = module_alpha (
#     pst = pst,
#     mod1 = mod1)
# res[[1]]
# res[[2]]


module_alpha = function(
    ps =ps,
    mod1 = mod1

){

  tem = mod1$group %>% table() %>%
    as.data.frame() %>%
    dplyr::arrange(desc(Freq))
  colnames(tem) = c("Model","OTU.num")
  head(tem)
  otu = NULL
  map = NULL
  for (i in 1:length(tem$Model)) {
    id.s = tem$Model %>% as.character()
    id.t =  mod1 %>% filter(group %in% id.s[i]) %>%.$ID
    ps.tem = subset_taxa.wt(ps %>% scale_micro(method = "sampling"), "OTU",id.t )

    otu = ps.tem %>% vegan_otu() %>% t() %>%
      as.data.frame()
    head(otu)
    colnames(otu) = paste(id.s[i],colnames(otu),sep = "_")
    map = data.frame(row.names = colnames(otu),ID = colnames(otu),Group = id.s[i])

    if (i == 1) {
      otu.f = otu
      map.f = map
    } else{

      otu$ID = row.names(otu)
      otu.f$ID = row.names(otu.f)
      tem.2 = otu.f %>% full_join(otu) %>%
        select(ID,everything()) %>%
        as.data.frame()
      row.names(tem.2) =  tem.2$ID
      tem.2$ID = NULL
      tem.2[is.na(tem.2)] = 0
      otu.f = tem.2
      map.f = rbind(map.f,map)
    }
  }

  pst.3 = phyloseq(
    otu_table(as.matrix(otu.f),taxa_are_rows = TRUE),
    sample_data(map.f) ,
    tax_table(ps)
  )
  index = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
  #--多种组合alpha分析和差异分析出图
  alp = ggClusterNet::alpha(ps = pst.3,inde="Shannon",group = "Group",Plot = TRUE,
              sampling = FALSE
  )
  index= alp
  head(index)

  #--提取三个代表指标作图
  sel = c(match("Shannon",colnames(index)),match("Richness",colnames(index)),match("Pielou_evenness",colnames(index)))
  data = cbind(data.frame(ID = 1:length(index$Group),group = index$Group),index[sel])
  head(data)

  # filename = paste(path,"Alpha_diversitydata",id3[m],".csv",sep = "")
  # write.csv(data,filename,quote = F)

  result = EasyStat::MuiaovMcomper2(data = data,num = c(3:5))

  # FileName <- paste(alppath,"/alpha_diversity_different_label.csv", sep = "")
  # write.csv(result,FileName,sep = "")
  # FileName <- paste(alppath,"/alpha_diversity_index.csv", sep = "")
  # write.csv(index,FileName,sep = "")
  sample_data(pst.3)
  result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = c(3:5),
                                            result = result,
                                            sig_show ="abc",ncol = 3 )
  p1_1 = result1[[1]] +
    ggplot2::guides(fill = guide_legend(title = NULL))

  p1_1

  plotdat = list(
    alpha = data,
    sigtab = result
  )


  #如何升级展示-提取数据用小提琴图展示
  p1_2 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
    geom_violin(alpha=1, aes(fill=group)) +
    geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
    labs(x="", y="")+
    facet_wrap(.~name,scales="free_y",ncol  = 3) +
    # theme_classic()+
    geom_text(aes(x=group , y=y ,label=stat)) +

    guides(color=guide_legend(title = NULL),
           shape=guide_legend(title = NULL),
           fill = guide_legend(title = NULL)
    )
  p1_2

  return(list(p1_1,p1_2,pst.3,plotdat))
}





# res = module_composition(pst = pst,mod1 = mod1,j = "Family")
# p1 = res[[1]]
# p2 = res[[2]]

module_composition = function(
    pst = pst,
    mod1 = mod1,
    j = "Family"
){
  tem = mod1$group %>% table() %>%
    as.data.frame() %>%
    dplyr::arrange(desc(Freq))
  colnames(tem) = c("Model","OTU.num")
  head(tem)

  otu = NULL
  map = NULL
  for (i in 1:length(tem$Model)) {
    id.s = tem$Model %>% as.character()
    id.t =  mod1 %>% filter(group %in% id.s[i]) %>%.$ID
    ps.tem = subset_taxa.wt(pst, "OTU", id.t )

    otu = ps.tem %>% vegan_otu() %>% t() %>%
      as.data.frame()
    head(otu)
    colnames(otu) = paste(id.s[i],colnames(otu),sep = "_")
    map = data.frame(row.names = colnames(otu),ID = colnames(otu),Group = id.s[i])

    if (i == 1) {
      otu.f = otu
      map.f = map
    } else{

      otu$ID = row.names(otu)
      otu.f$ID = row.names(otu.f)
      tem.2 = otu.f %>% full_join(otu) %>%
        select(ID,everything()) %>%
        as.data.frame()
      row.names(tem.2) =  tem.2$ID
      tem.2$ID = NULL
      tem.2[is.na(tem.2)] = 0
      otu.f = tem.2
      map.f = rbind(map.f,map)
    }
  }

  pst.2 = phyloseq(
    otu_table(as.matrix(otu.f),taxa_are_rows = TRUE),
    sample_data(map.f) ,
    tax_table(pst)
  )
  # library(RColorBrewer)
  # colset1 <- brewer.pal(10,"Paired")
  # pst.2 = pst.2  %>%
  #   subset_taxa(
  #     !Genus %in% "Unassigned"
  #   )

  result = barMainplot(ps = pst.2,
                       tran = F,
                       j = j,
                       # axis_ord = axis_order,
                       label = FALSE,
                       sd = FALSE,
                       Top = 10)
  p4_1 <- result[[1]] +
    scale_fill_hue() + theme_classic()
  p4_1
  p4_2  <- result[[3]]  +
    scale_fill_hue() + theme_classic()
  p4_2

  tem1 = result[[2]]

  result = barMainplot(ps = pst.2,
                       tran = T,
                       j = j,
                       # axis_ord = axis_order,
                       label = FALSE,
                       sd = FALSE,
                       Top = 10)
  p3_1 <- result[[1]] +
    scale_fill_hue() + theme_classic()
  p3_1
  p3_2  <- result[[3]]  +
    scale_fill_hue() + theme_classic()
  p3_2
  tem2 = result[[2]]
  library(patchwork)
  p00 = p4_1|p3_1
  p01 = p4_2|p3_2

  plotdat = list(
    bundance = tem1,
    relaabundance = tem2

  )

  return(list(p00,p01,pst.2,plotdat = plotdat))
}







# res  = module_display(
#     pst = pst,
#     r.threshold= 0.6,
#     p.threshold=0.1,
#     select.mod = c("model_1","model_2","model_3","model_4"))
# res[[1]]
# res[[2]]
module_display = function(
    pst = pst,
    Top = 500,
    r.threshold= 0.8,
    p.threshold=0.05,
    select.mod = c("model_1","model_2","model_3")
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
  result = model_maptree2(cor = cor, method = "cluster_walktrap" )
  node = result[[1]]
  netClu = result[[2]]
  branch = result[[3]] %>% filter(!str_detect(elements,"other"))
  head(branch)
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
  # head(bra)

  #---选择模块的微生物网络
  # select.mod = 3
  # select.mod = select.mod
  # mod1 = result[[2]]
  # result = model_maptree2(cor = cor, method = "cluster_walktrap" )
  head(mod1)
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
    # geom_text(aes(x = x, y = y,label = elements), data = bra) +
    # labs( title = paste(layout,"network",sep = "_"))+
    # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
    # discard default grid + titles in ggplot2
    theme(panel.background = element_blank()) +
    # theme(legend.position = "none") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  # pnet


  return(list(pnet,p2,cor,netClu))
}
