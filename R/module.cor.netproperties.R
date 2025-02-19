
module.cor.netproperties = function(
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
  dat.f = netproperties.sample(pst = pst,cor = cor)
  # head(dat.f)
  dat.f$id = row.names(dat.f)
  dat.f = dat.f %>% dplyr:: select(id,everything())

  colnames(env)[1] = "id"
  subenv = env %>% dplyr::select(id,everything()) %>% dplyr::select(id,select.env )
  # head(data)

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










  plotdat = list(
    preopertites.env = preoptab

  )
  return(list(p0_1,dat.f,plotdat))

}

