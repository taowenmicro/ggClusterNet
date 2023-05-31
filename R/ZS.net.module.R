



# dat = ZS.net.module(
#     pst = pst,
#     Top = 500,
#     r.threshold= 0.8,
#     p.threshold=0.05,
#     select.mod =  NULL
#     )
#
# head(dat)




ZS.net.module = function(
    pst = pst,
    Top = 500,
    corg = NULL,
    method = "spearman",
    r.threshold= 0.8,
    p.threshold=0.05,
    select.mod = c("model_1","model_2","model_3")

){

  # result = cor_Big_micro(ps = pst,
  #                        N = Top,
  #                        r.threshold= r.threshold,
  #                        p.threshold= p.threshold,
  #                        method = "spearman")
  # cor = result[[1]]

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

  } else if (is.character(select.mod)) {
    select.mod.name = select.mod
    mod1 = mod1 %>% filter(!group == "mother_no",
                           group %in%c(select.mod.name)

    ) %>% select(ID,group,degree)

  } else if(is.null(select.mod)){
    mod1 = mod1
  }

  #-计算每个模块的平均Zscore
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


  return(Zscore = tem)
}
