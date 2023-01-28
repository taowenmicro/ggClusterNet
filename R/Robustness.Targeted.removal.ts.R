

# res = Robustness.Targeted.removal.ts(
#     ps.st,
#     N = 200,
#     degree = TRUE,
#     zipi = FALSE,
#     r.threshold= 0.8,
#     p.threshold=0.05,
#     method = "spearman",
#     order = "time",
#     g1 = "Group",# 分组1
#     g2 = "space",# 分组2
#     g3 = "time")


Robustness.Targeted.removal.ts = function(
    ps.st,
    N = 200,
    degree = TRUE,
    zipi = FALSE,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman",
    order = "time",
    g1 = "Group",# 分组1
    g2 = "space",# 分组2
    g3 = "time",# 分组3
    ord.g1 =NULL,# 排序顺序
    ord.g2 = NULL, # 排序顺序
    ord.g3 = NULL# 排序顺序
    ){

  otutab<- ps.st %>%
    vegan_otu() %>%
    as.data.frame()
  dim(otutab)

  # id <- sample_data(ps)$Group %>% unique()
  # i  = 1
  #计算每个物种的平均丰度，使用测序深度标准化
  sp.ra<-colMeans(otutab)/mean(rowSums(otutab))   #relative abundance of each species


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
      temx = paste(dat.f[j,1],dat.f[j,2],group[n],sep = ".")
      cor.all[[temx]] = cor
    #存在某些情况计算不出来相关系数，定义相关为0
    cor[is.na(cor)]<-0
    #-去除自相关的点
    diag(cor)<-0
    #-查看网络边的数量
    sum(abs(cor)>0)/2
    #网络中节点的数量
    sum(colSums(abs(cor))>0)  # node number: number of species with at least one linkage with others.
    #去除没有任何相关的节点.
    network.raw<-cor[colSums(abs(cor))>0,colSums(abs(cor))>0]

    ##read otu table

    otutab<- ps.g %>%
      scale_micro() %>%
      subset_taxa.wt(
        "OTU", row.names(network.raw)
      ) %>%
      vegan_otu() %>%
      t() %>%
      as.data.frame()

    #对应的删除otu表格otu
    sp.ra2<- rowSums(otutab)
    # sp.ra2
    sum(row.names(network.raw) %in% names(sp.ra2))  #check if matched
    ## robustness simulation
    igraph = make_igraph(cor)

    module.hub = NULL
    if (zipi ) {
      res = ZiPiPlot(igraph = igraph,method = "cluster_fast_greedy")
      # p <- res[[1]]
      model = res[[2]] %>% filter(roles == "Module hubs")
      head(model)
      model$roles %>% unique()
      module.hub <- as.character(row.names(model))
    }

    if (length( module.hub) == 0) {
      degree = TRUE
    }

    if (degree) {
      ret3 = node_properties(igraph) %>%
        as.data.frame() %>%
        filter(!is.na(igraph.degree) ) %>%
        arrange(desc(igraph.degree))
      # head(ret3)
      tem = round(length(ret3$igraph.degree) * 0.05,0)
      module.hub = row.names(ret3)[1:tem]
    }


    rmsimu<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){

      t(sapply(rm.p.list,function(x){
        remains=sapply(1:nperm,function(i){
          rand.remov2.once(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
        })
        remain.mean=mean(remains)
        remain.sd=sd(remains)
        remain.se=sd(remains)/(nperm^0.5)
        result<-c(remain.mean,remain.sd,remain.se)
        names(result)<-c("remain.mean","remain.sd","remain.se")
        result
      }))
    }


    Weighted.simu<-rmsimu(netRaw=network.raw,
                          rm.p.list=1:length(module.hub),
                          keystonelist=module.hub,
                          sp.ra=sp.ra2,
                          abundance.weighted=T,nperm=100)
    Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub),
                            keystonelist=module.hub,
                            sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

    dat1<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),
                     rbind(Weighted.simu,Unweighted.simu),
                     weighted=rep(c("weighted","unweighted"),
                                  each=length(module.hub)),
                     Group= temx)

    if (j ==1 & n == 1) {
      dat.t = dat1
    } else{
      dat.t = rbind(dat.t,dat1)
    }

    }
  }


  #-行--空间
  if (order == "space"|order == "g2") {
    row.id = g3
    row.num =  length(ti)
    a = c()
    for (i in 1:length(sp)) {
      for (j in 1:length(ti)) {
        tem = paste(sp[i],ti[j],group,sep = ".")
        a = c(a,tem)
      }
    }

  } else if (order == "time"|order == "g3") {
    row.id = g2
    row.num =  length(sp)
    a = c()
    for (j in 1:length(ti)) {
      for (i in 1:length(sp)) {
        tem = paste(sp[i],ti[j],group,sep = ".")
        a = c(a,tem)
      }
    }
  }

  head(dat.t)
  dat.t$group = sapply(strsplit(as.character(dat.t$Group), "[.]"), `[`, 3)
  dat.t$time = sapply(strsplit(as.character(dat.t$Group), "[.]"), `[`, 2)
  dat.t$space = sapply(strsplit(as.character(dat.t$Group), "[.]"), `[`, 1)
  dat.t$label = paste(dat.t$time,dat.t$space,sep = ".")
  dat.t$Group = factor(dat.t$Group,levels = a)

  p = ggplot(dat.t[dat.t$weighted=="weighted",],
             aes(x=Number.hub.removed, y=remain.mean, group=group, color=group)) +
    geom_line()+
    geom_pointrange(aes(ymin=remain.mean-remain.sd,
                        ymax=remain.mean+remain.sd),size=0.2)+
    facet_wrap(.~ label,scales="free_y",ncol = row.num ) +
    xlab("Number of module hubs removed")+
    ylab("Proportion of species remained")+
    theme_light()


  p2 = ggplot(dat.t[dat.t$weighted=="unweighted",],
              aes(x=Number.hub.removed, y=remain.mean, group=group, color=group)) +
    geom_line()+
    geom_pointrange(aes(ymin=remain.mean-remain.sd,
                        ymax=remain.mean+remain.sd),size=0.2)+
    facet_wrap(.~ label,scales="free_y",ncol = row.num ) +
    xlab("Number of module hubs removed")+
    ylab("Proportion of species remained")+
    theme_light()

  plot = list(weighted = p,
              unweighted = p2
  )

  return(list(plot,dat.t))
}




