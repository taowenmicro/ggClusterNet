# library(ggClusterNet)
# library(phyloseq)
# library(tidyverse)


# dat = module.compare.net.pip(
#     ps = ps,Top = 500,degree = TRUE,
#     zipi = FALSE,r.threshold= 0.8,
#     p.threshold=0.05,
#     method = "spearman",padj = F,n = 3)
# res = dat[[1]]
# head(res)
#
# #--cor matrix相关矩阵提取
# cor = dat[[2]]
# names(cor)
#
# #--分组展示网络中的节点
# edgtb = dat[[3]]
# head(edgtb)
# edgtb$group %>% table()
#


module.compare.net.pip = function(
    ps = ps,
    Top = NULL,
    corg = NULL,
    degree = TRUE,
    zipi = FALSE,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman",
    padj = F,
    n = 3
){

  # map = sample_data(ps)
  # head(map)
  # id <- map$Group %>% unique()
  # otu = ps %>% vegan_otu() %>% t() %>%
  #   as.data.frame()
  # head(otu)
  # tax = ps %>% vegan_tax() %>%
  #   as.data.frame()
  # head(tax)
  if (!is.null(corg)) {
    id = names(corg)
  } else if (is.null(corg)){

  }

  cortb = list()
  for (i in 1:length(id)) {

    if (is.null(corg)) {
      pst =  ps %>%
        scale_micro() %>%
        subset_samples.wt("Group", c(id[i])) %>%
        filter_OTU_ps(Top)

      result = cor_Big_micro(ps = pst,
                             N = 0,
                             r.threshold= r.threshold,
                             p.threshold= p.threshold,
                             method = method)

      cor = result[[1]]
      # head(cor)
    } else if (!is.null(corg)){
      cor = corg[[id[i]]]
    }


    cortb[[i]] = cor
    names(cortb)[[i]] = id[i]
    # igraph = make_igraph(cor)
    #--计算模块信息，部分OTU没有模块，注意去除
    result2 = model_maptree2(cor = cor,
                             method = "cluster_fast_greedy"
    )

    mod1 = result2[[2]]
    head(mod1)
    # mod1$group =  NULL
    mod1 = mod1 %>%
      dplyr::filter(!group == "mother_no") %>%
      dplyr::select(ID,group)
    # mod1$group = paste(id[i],mod1$group,sep = "")
    mod1$Group = id[i]
    mod1$group = mod1$Group
    head(mod1)

    if (i == 1) {
      dat = mod1
    } else {
      dat = rbind(dat,mod1)
    }
  }
  # library(tidyfst)
  node_table2  = dat
  head(node_table2)
  node_table2$Group %>% table()
  head(dat)
  dat2 = model_compare.net(
    node_table2 = dat,
    n = n,
    padj = TRUE
  )

  return(list(dat2,
              cortab =  cortb,
              edgetab = dat

              ))
}


model_compare.net = function(
    node_table2 = node_table2,
    n = 3,
    padj = F
){
  node_table2$value = 1
  #提取全部的模块及其名称
  module_list = unique(node_table2$group)
  map = node_table2[,2:3] %>% distinct(group, .keep_all = TRUE)
  mytable = node_table2 %>% tidyfst::df_mat(ID,group,value) %>% as.matrix()
  mytable[is.na(mytable)] <- 0
  #--去除少于n个OTU的模块
  if (colSums(mytable)[colSums(mytable)< n] %>% length() == 0) {
    mytable_kp = mytable
  } else {
    mytable_kp = mytable[,-which(colSums(mytable)<n)]
  }

  # 对应的处理一下map文件
  head(map)
  map_kp = map %>% filter(group %in% colnames(mytable_kp))
  mytable_kp = as.data.frame(mytable_kp)

  #---构造分组两两组合
  id = unique(node_table2$Group)

  network_pair = combn(id,2) %>% t() %>%
    as.matrix()


  total_mod_pairs = matrix(NA, nrow=nrow(network_pair), ncol=3)
  # i = 1
  for (i in 1:nrow(network_pair)){
    # 对两个组的每一个模块都进行比对
    module_pair = as.matrix(expand.grid(
      map_kp$group[which(map_kp$Group==network_pair[i,1])],
      map_kp$group[which(map_kp $Group==network_pair[i,2])]))

    total_mod_pairs[i,] = c(network_pair[i,], nrow(module_pair))
  }

  sig_mod_pairs = matrix(NA, nrow=0, ncol=4)
  sig_detailed_table = c("module1", "module2", "both", "P1A2", "P2A1", "A1A2", "p_raw", "p_adj")

  for (i in 1:nrow(network_pair)){
    # 全部的需要比对的模块
    module_pair = as.matrix(expand.grid(map_kp$group[which(map_kp$Group==network_pair[i,1])],
                                        map_kp$group[which(map_kp $Group==network_pair[i,2])]))
    overlap = apply(module_pair, 1, FUN= find_overlap, bigtable= mytable_kp)
    only1 = apply(module_pair, 1, FUN= find_only_in_1, bigtable= mytable_kp)
    only2 = apply(module_pair, 1, FUN= find_only_in_2, bigtable= mytable_kp)
    denominator = apply(module_pair, 1, FUN= find_N2, mapping=map, bigtable= mytable)
    none = denominator-(overlap + only1 + only2)
    count_table = data.frame(module1 = module_pair[,1],
                             module2 = module_pair[,2],
                             Both=overlap,
                             P1A2=only1,
                             P2A1=only2, A1A2=none)

    p_raw=c()
    # tt = 1
    for (tt in 1:nrow(count_table))
    {
      x=count_table[tt,]
      p = fisher_test(x)
      p_raw = c(p_raw, p)
    }

    count_table$p_raw = p_raw
    if (padj){
      count_table$p_adj = p.adjust(count_table$p_raw, method = "bonferroni")
    }else{
      count_table$p_adj = count_table$p_raw
    }

    if (i == 1) {
      dat = count_table
    } else {
      dat = rbind(dat,count_table)
    }
  }


  return(dat)
}

find_N2 = function(mods, mapping, bigtable){
   match_nwk1_nwk2 = dim(bigtable)[1]
  return(match_nwk1_nwk2)
}

