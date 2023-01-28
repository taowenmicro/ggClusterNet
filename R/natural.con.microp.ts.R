

# res = natural.con.microp.ts(
#   ps.st,
#   N = 200,
#   r.threshold= 0.8,
#   p.threshold=0.05,
#   method = "spearman",
#   order = "time",
#   g1 = "Group",# 分组1
#   g2 = "space",# 分组2
#   g3 = "time",# 分组3
#   ord.g1 =NULL,# 排序顺序
#   ord.g2 = NULL, # 排序顺序
#   ord.g3 = NULL,# 排序顺序
#   norm = F,
#   end = N -2,
#   start = 0,
#   con.method = "pulsar"
#
#   )
#
# res[[1]]

natural.con.microp.ts = function(
    ps.st,
    N = 200,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman",
    order = "time",
    g1 = "Group",# 分组1
    g2 = "space",# 分组2
    g3 = "time",# 分组3
    ord.g1 =NULL,# 排序顺序
    ord.g2 = NULL, # 排序顺序
    ord.g3 = NULL,# 排序顺序
    norm = F,
    end = N -2,
    start = 0,
    con.method = "pulsar"

){

  otutab<- ps.st %>%
    vegan_otu() %>%
    as.data.frame()
  dim(otutab)

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




    dat = natural.con.micro(cor = cor,
                            start = start,
                            end = 198,
                            norm = norm,
                            method = con.method
    )
    dat$Group = as.character(temx)
    if (j ==1 & n == 1) {
      dat.t = dat
    }else{
      dat.t = rbind(dat.t,dat)
    }

    }
  }
  names(cor.all)
  head(dat.t)



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
  dat.t$Group %>% unique()
  head(dat.t)
  dat.t$group = sapply(strsplit(as.character(dat.t$Group), "[.]"), `[`, 3)
  dat.t$time = sapply(strsplit(as.character(dat.t$Group), "[.]"), `[`, 2)
  dat.t$space = sapply(strsplit(as.character(dat.t$Group), "[.]"), `[`, 1)
  dat.t$label = paste(dat.t$time,dat.t$space,sep = ".")
  dat.t$Group = factor(dat.t$Group,levels = a)

  p=ggplot(dat.t, aes(Num.of.remove.nodes, Natural.connectivity,color=group)) +
    geom_point(alpha = 0.1) +
    geom_smooth(level = F) +
    facet_wrap(.~ label,scales="free_y",ncol = row.num ) +
    theme_set(theme_bw(base_size = 25))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black")) + theme_bw()
  p
  # p
  return(list(p,dat.t))
}
