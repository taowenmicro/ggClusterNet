
# res = community.stability.ts (
#     ps.st = ps.st,
#     N = 200,
#     r.threshold= 0.8,
#     p.threshold=0.05,
#     method = "spearman",
#     order = "time",
#     g1 = "Group",# 分组1
#     g2 = "space",# 分组2
#     g3 = "time",# 分组3
#     map.art = NULL, # 人工输入的分组 默认为NULL
#     time = F,# 稳定性是否有时间序列
#     ord.map = TRUE# map文件是否是已经按照pair要求进行了排序
# )
#
# res[[1]]
# res[[2]]


#
community.stability.ts = function(
  ps.st = ps.st,
  N = 200,
  r.threshold= 0.8,
  p.threshold=0.05,
  method = "spearman",
  order = "space",
  g1 = "Group",# 分组1
  g2 = "space",# 分组2
  g3 = "time",# 分组3
  ord.g1 =NULL,# 排序顺序
  ord.g2 = NULL,# 排序顺序
  ord.g3 = NULL,# 排序顺序
  map.art = NULL, # 人工输入的分组 默认为NULL
  time = F,
  ord.map = TRUE# map文件是否是已经按照pair要求进行了排序
  ){


  otutab<- ps.st %>%
    vegan_otu() %>%
    t() %>%
    as.data.frame()
  dim(otutab)

  ps.all = ps.st
  map = sample_data(ps.all)

  treat = ps.st %>% sample_data()
  treat$ID = row.names(treat)
  head(treat)

  tem1 = treat[,g1] %>% as.matrix() %>% as.vector()
  tem2 = treat[,g2] %>% as.matrix() %>% as.vector()
  tem3 = treat[,g3] %>% as.matrix() %>% as.vector()

  tem4 = paste(tem3,tem2,tem1,sep = ".")
  tem5 = tem4 %>% table() %>%as.data.frame() %>% .$Freq %>% unique() %>% length()
  rep = tem4 %>% table() %>%as.data.frame() %>% .$Freq %>% unique()
  num.g = unique(tem4) %>% length()
  #-如果重复数量相同，并且重复排序相同
  if (tem5 == 1) {
    treat$allg = tem4
    head(treat)
    if (ord.map == F) {
      treat = treat %>%as.tibble() %>%arrange(desc(allg)) %>% as.data.frame()

    }else if(ord.map == T) {#---重复数量相同，使用原来顺序制作pair
      treat = treat
    }
    treat$pair = paste( "A",c(rep(1:rep,num.g)),sep = "")
    row.names(treat) = treat$ID
    # sample_data(ps.st) = treat
  } else if (tem5 > 1) {
    #--如果重复数量不等，或者部分相同部分不同
    # 选择可利用的进行分析-这部分暂时不进行书写，分析者提供子map文件进行分析
    print("Repeats not the same number")
    treat = map.art
  }
  sample_data(ps.st) = treat
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


  comm = otutab  %>% t()
  #去除NA值
  sum(is.na(comm)) # check NA
  comm[is.na(comm)]=0# if have, should change to zero
  plot.lev = unique(treat$pair)

  # time = F
  #-提取时间序列
  year.lev = sort(unique(treat$allg))

  if (time == TRUE) {
    #-构造序列
    zeta.lev = 2:length(year.lev)

    # 构造从2到6的全部这组合，这里使用断棍模型构造全部组合
    year.windows=lapply(1:length(zeta.lev),
                        function(i)
                        {zetai=zeta.lev[i]
                        lapply(1:(length(year.lev)-zetai+1),function(j){year.lev[j:(j+zetai-1)]})
                        })

    names(year.windows)=zeta.lev
    year.windows
  } else if(time == FALSE){

    tem2 = list()
    A = c()
    tem = combn(year.lev ,2)
    # tem
    for (i in 1:dim(tem)[2]) {
      tem2[[i]] = c(tem[1,i],tem[2,i])
      A[i]= paste("Zeta",tem[1,i],tem[2,i],sep = "_")
    }
    names(tem2) = A
    zeta.lev = rep(2,length(A))
    year.windows = list()
    year.windows[[1]] = tem2
    names(year.windows) = "2"
  }



  # year.windows %>% names() %>% length()


  # 基于不同分组样本的群落稳定性功能函数:物种最小丰度和乘以样本数量，
  # 得到的结果除以多组全部微生物丰度的和

  comstab<-function(subcom){((nrow(subcom)*sum(apply(subcom,2,min)))/sum(subcom))^0.5}
  # subcom = comijk
  i = 1
  j = 1
  k = 1
   stabi=t(sapply(1:length(plot.lev),
                                function(j)
                                {
                                  plotj=plot.lev[j]
                                  sapply(1:length(year.windows[[i]]),
                                         function(k)
                                         {
                                           yearwdk=year.windows[[i]][[k]] %>% as.character()
                                           sampijk=rownames(treat)[which((treat$pair==plotj) & (treat$allg %in% yearwdk))]
                                           outijk=NA
                                           if(length(sampijk) < length(yearwdk))
                                           {
                                             warning("plot ",plotj," has missing year in year window ",paste(yearwdk,collapse = ","))
                                           }else if(length(sampijk) > length(yearwdk)){
                                             warning("plot ",plotj," has duplicate samples in at least one year of window ",paste(yearwdk,collapse = ","))
                                           }else{
                                             comijk=comm[which(rownames(comm) %in% sampijk),,drop=FALSE]
                                             outijk=comstab(comijk)
                                           }
                                           outijk
                                         })
                                }))
                 if(nrow(stabi)!=length(plot.lev) & nrow(stabi)==1){stabi=t(stabi)}
                 rownames(stabi) = plot.lev
                 colnames(stabi)=sapply(year.windows[[i]],function(v){paste0("Zeta",zeta.lev[i],paste0(v,collapse = "_"))})
                 # stabi
               # })


  head(stabi)
  stabl = list()
  stabl[[1]] = stabi

  stabm=Reduce(cbind,stabl) %>% as.data.frame()
  head(stabm)
  dat = stabm %>%
    # rownames_to_column("id") %>%
    gather()
  head(dat)
  dim(dat)[1]/6

  #-行--空间
  if (order == "space"|order == "g2") {
    row.id = g3
    row.num =  length(ti)
    a = c()
    for (i in 1:length(sp)) {
      for (j in 1:length(ti)) {
        tem = paste(sp[i],ti[j],sep = ".")
        a = c(a,tem)
      }
    }

  } else if (order == "time"|order == "g3") {
    row.id = g2
    row.num =  length(sp)
    a = c()
    for (j in 1:length(ti)) {
      for (i in 1:length(sp)) {
        tem = paste(sp[i],ti[j],sep = ".")
        a = c(a,tem)
      }
    }
  }


  dat$vs1= sapply(strsplit(as.character(dat$key), "[_]"), `[`, 1)
  dat$vs1 = gsub("Zeta2","",dat$vs1)
  dat$vs1.t= sapply(strsplit(as.character(dat$vs1), "[.]"), `[`, 1)
  dat$vs1.s= sapply(strsplit(as.character(dat$vs1), "[.]"), `[`, 2)
  dat$vs1.g= sapply(strsplit(as.character(dat$vs1), "[.]"), `[`, 3)

  dat$vs2 = sapply(strsplit(as.character(dat$key), "[_]"), `[`, 2)
  dat$vs2.t= sapply(strsplit(as.character(dat$vs2), "[.]"), `[`, 1)
  dat$vs2.s= sapply(strsplit(as.character(dat$vs2), "[.]"), `[`, 2)
  dat$vs2.g= sapply(strsplit(as.character(dat$vs2), "[.]"), `[`, 3)

  dat$ts1 = paste(dat$vs1.t,dat$vs1.s,sep = "_")
  dat$ts2 = paste(dat$vs2.t,dat$vs2.s,sep = "_")
  dat$group = paste(dat$vs1.g,dat$vs2.g,sep = "_")

  # dat2$crosstype  %>% unique()
  # head(dat)
  dat2 = dat %>% dplyr::mutate(
                               crosstype = ifelse(ts1 == ts2, as.character(ts2), "across")) %>%
    filter(crosstype != "across")
  # head(dat2)
  # dat.t$group = sapply(strsplit(as.character(dat.t$), "[.]"), `[`, 3)
  dat2$time = sapply(strsplit(as.character(dat2$crosstype), "[_]"), `[`, 1)
  dat2$space = sapply(strsplit(as.character(dat2$crosstype), "[_]"), `[`, 2)
  dat2$label = paste(dat2$space,dat2$time,sep = ".")
  dat2$label = factor(dat2$label,levels = a)
  # head(dat2)


  p = ggplot(dat2) + geom_boxplot(aes(x = group,y = value,fill = group)) + theme_bw() +
    labs(y = "community.stability") +
    facet_wrap(.~ label,scales="free_y",ncol = row.num )
  p



  return(list(p,dat2))
}

