
#--模块比较

# res = module.compare.m(
#   ps = ps,
#   Top = 500,
#   degree = TRUE,
#   zipi = FALSE,
#   r.threshold= 0.8,
#   p.threshold=0.05,
#   method = "spearman",
#   padj = F,
#   n = 3)
# p = res[[1]]
# dat = res[[2]]
# head(dat)
# dat2 = res[[3]]
# head(dat2)


module.compare.m = function(
    ps = ps,
    Top = 500,
    degree = TRUE,
    zipi = FALSE,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman",
    padj = F,
    n = 3

){
  map = sample_data(ps)
  head(map)
  id <- map$Group %>% unique()
  otu = ps %>% vegan_otu() %>% t() %>%
    as.data.frame()
  head(otu)
  tax = ps %>% vegan_tax() %>%
    as.data.frame()
  head(tax)


  #-判断相同模块#----
  i= 1

  for (i in 1:length(id)) {

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
    head(cor)

    # igraph = make_igraph(cor)
    #--计算模块信息，部分OTU没有模块，注意去除
    result2 = model_maptree2(cor = cor,
                             method = "cluster_fast_greedy"
    )

    mod1 = result2[[2]]
    head(mod1)

    mod1 = mod1 %>% filter(!group == "mother_no") %>% select(ID,group)
    mod1$group = paste(id[i],mod1$group,sep = "")
    mod1$Group = id[i]
    head(mod1)

    if (i == 1) {
      dat = mod1
    } else {
      dat = rbind(dat,mod1)
    }
  }

  node_table2  = dat
  head(node_table2)
  node_table2$Group %>% table()
  library(tidyfst)
  dat2 = model_compare(
    node_table2 = dat,
    n = n,
    padj = padj
  )

  head(dat2)
  head(node_table2)
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
                                mi.size = 0.5,
                                zoom = 0.2)
  node = result2[[1]]
  head(node)
  head(edge)

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
    ggrepel::geom_text_repel(aes(X1, X2,label=elements),size=4, data = node)+
    # discard default grid + titles in ggplot2
    theme(panel.background = element_blank()) +
    # theme(legend.position = "none") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  pnet
  return(list(pnet,node_table2,dat2))
}

# res= Robustness.Targeted.removal(ps = ps,
#                                  Top = 500,
#                                  degree = TRUE,
#                                  zipi = FALSE,
#                                  r.threshold= 0.8,
#                                  p.threshold=0.05,
#                                  method = "spearman")
#
# p = res[[1]]
# p
# dat = res[[2]]
Robustness.Targeted.removal = function(
    ps = ps,
    Top = 500,
    degree = TRUE,
    zipi = FALSE,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman"

){

  otutab<- ps %>%
    vegan_otu() %>%
    as.data.frame()
  dim(otutab)

  id <- sample_data(ps)$Group %>% unique()
  i  = 1
  #计算每个物种的平均丰度，使用测序深度标准化
  sp.ra<-colMeans(otutab)/mean(rowSums(otutab))   #relative abundance of each species

  for (i in 1:length(id)){
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
    head(cor)


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

    otutab<- pst %>%
      scale_micro() %>%
      subset_taxa.wt(
        "OTU", row.names(network.raw)
      ) %>%
      vegan_otu() %>%
      t() %>%
      as.data.frame()



    #对应的删除otu表格otu
    sp.ra2<- rowSums(otutab)
    sp.ra2
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
      head(ret3)
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
                     Group= id[i])

    if (i ==1) {
      dat.f = dat1
    } else if(i != 1){
      dat.f = rbind(dat.f,dat1)
    }

  }


  p = ggplot(dat.f[dat.f$weighted=="weighted",],
             aes(x=Number.hub.removed, y=remain.mean, group=Group, color=Group)) +
    geom_line()+
    geom_pointrange(aes(ymin=remain.mean-remain.sd,
                        ymax=remain.mean+remain.sd),size=0.2)+

    xlab("Number of module hubs removed")+
    ylab("Proportion of species remained")+
    theme_light()


  p2 = ggplot(dat.f[dat.f$weighted=="unweighted",],
              aes(x=Number.hub.removed, y=remain.mean, group=Group, color=Group)) +
    geom_line()+
    geom_pointrange(aes(ymin=remain.mean-remain.sd,
                        ymax=remain.mean+remain.sd),size=0.2)+

    xlab("Number of module hubs removed")+
    ylab("Proportion of species remained")+
    theme_light()

  p3 = p | p2

  return(list(p3,dat.f,p,p2))
}



#---网络易损性
# res = Vulnerability.micro(ps = ps,
#                           Top = 500,
#                           degree = TRUE,
#                           zipi = FALSE,
#                           r.threshold= 0.8,
#                           p.threshold=0.05,
#                           method = "spearman")
#
#
# p = res[[1]]
# dat = res[[2]]
Vulnerability.micro = function(
    ps = ps,
    Top = 500,
    degree = TRUE,
    zipi = FALSE,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman"
){

  otutab<- ps %>%
    vegan_otu() %>%
    as.data.frame()
  dim(otutab)

  id <- sample_data(ps)$Group %>% unique()
  i  = 1
  A = c()

  for (i in 1:length(id)){
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
    head(cor)

    vulnerability = function(cor = cor){

      cor[abs(cor)>0]<-1 # adjacency matrix
      g = graph_from_adjacency_matrix(as.matrix(cor),
                                      mode="undirected",
                                      weighted = NULL, diag = FALSE,
                                      add.colnames = NULL) # note: this graph contains isolated nodes.
      # remove isolated nodes
      iso_node_id = which(igraph::degree(g)==0)
      g2 = igraph::delete.vertices(g, iso_node_id) # graph without isolated nodes

      #check node number and links
      length(igraph::V(g2));length(igraph::E(g2))

      # calculate vulnerability of each node
      node.vul<-info.centrality.vertex(g2)
      return(max(node.vul))
    }
    A[i] = vulnerability(cor = cor)
  }
  dat = data.frame(ID = id,Vulnerability = A)
  head(dat)
  p = ggplot(dat) + geom_bar(aes(x = ID,y = Vulnerability,fill = ID),stat =
                               "identity")
  return(list(p,dat))
}



#--计算负相关的比例

# res = negative.correlation.ratio(ps = ps,
#                                  Top = 500,
#                                  degree = TRUE,
#                                  zipi = FALSE,
#                                  r.threshold= 0.8,
#                                  p.threshold=0.05,
#                                  method = "spearman")
#
#
# p = res[[1]]
# dat = res[[2]]
negative.correlation.ratio = function(
    ps = ps,
    Top = 500,
    degree = TRUE,
    zipi = FALSE,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman"
){

  otutab<- ps %>%
    vegan_otu() %>%
    as.data.frame()
  dim(otutab)

  id <- sample_data(ps)$Group %>% unique()
  i  = 1
  B = c()
  for (i in 1:length(id)) {
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
    head(cor)
    igraph = make_igraph(cor)
    ret2 = net_properties.2(igraph) %>% as.data.frame()
    head(ret2)
    a = ret2[1,1] %>% as.numeric()
    n = ret2[3,1] %>% as.numeric()
    B[i] = n/a *100
  }


  dat = data.frame(ID = id,ratio = B)

  p = ggplot(dat) + geom_bar(aes(x = ID,y = ratio,fill = ID),stat =
                               "identity")
  return(list(p,dat))
}

#-群落稳定性
# treat = ps %>% sample_data()
# treat$pair = paste( "A",c(rep(1:6,3)),sep = "")
# head(treat)
# sample_data(ps) = treat
# res = community.stability( ps = ps)
#
# p = res[[1]]
# dat = res[[2]]
community.stability = function(
    ps = ps,
    time = TRUE
){
  # id <- sample_data(ps)$Group %>% unique()
  otutab = ps %>%
    scale_micro() %>%
    # phyloseq::subset_samples(Group %in% c(id[i])) %>%
    vegan_otu() %>% t() %>%
    as.data.frame()
  head(otutab)
  comm=otutab %>% t()
  #去除NA值
  sum(is.na(comm)) # check NA
  comm[is.na(comm)]=0# if have, should change to zero



  plot.lev=unique(treat$pair)
  #-提取时间序列
  year.lev = sort(unique(treat$Group))
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
    for (i in 1:length(year.lev)) {
      tem2[[i]] = c(tem[1,i],tem[2,i])
      A[i]= paste("Zeta",tem[1,i],tem[2,i],sep = "_")
    }
    names(tem2) = A

    zeta.lev = rep(2,length(A))

  }

  year.windows = list()
  year.windows[[1]] = tem2
  names(year.windows) = "2"

  year.windows


  # 基于不同分组样本的群落稳定性功能函数物种最小丰度和乘以样本数量，得到的结果除以多组全部微生物丰度的和
  comstab<-function(subcom){((nrow(subcom)*sum(apply(subcom,2,min)))/sum(subcom))^0.5}
  # subcom = comijk

  stabl=lapply(1:length(year.windows),
               function(i)
               {
                 stabi=t(sapply(1:length(plot.lev),
                                function(j)
                                {
                                  plotj=plot.lev[j]
                                  sapply(1:length(year.windows[[i]]),
                                         function(k)
                                         {
                                           yearwdk=year.windows[[i]][[k]] %>% as.character()
                                           sampijk=rownames(treat)[which((treat$pair==plotj) & (treat$Group %in% yearwdk))]
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
                 colnames(stabi)=sapply(year.windows[[i]],function(v){paste0("Zeta",zeta.lev[i],paste0(v,collapse = ""))})
                 stabi
               })




  stabm=Reduce(cbind,stabl) %>% as.data.frame()
  dat = stabm %>%
    # rownames_to_column("id") %>%
    gather()
  p = ggplot(dat) + geom_boxplot(aes(x = key,y = value,fill = key)) + theme_bw() +
    labs(y = "community.stability")

  return(list(p,dat))
}


# res = Robustness.Random.removal(ps = ps,
#                                 Top = 500,
#                                 r.threshold= 0.8,
#                                 p.threshold=0.05,
#                                 method = "spearman")
# p = res[[1]]
# dat = res[[2]]
Robustness.Random.removal = function(
    ps = ps,
    Top = 500,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman"
){
  otutab<- ps %>%
    vegan_otu() %>%
    as.data.frame()
  dim(otutab)

  id <- sample_data(ps)$Group %>% unique()
  i  = 1
  #计算每个物种的平均丰度，使用测序深度标准化
  sp.ra<-colMeans(otutab)/mean(rowSums(otutab))   #relative abundance of each species

  # library(ggClusterNet)
  for (i in 1:length(id)){
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
    head(cor)

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
    #对应的删除otu表格otu
    sp.ra2<-sp.ra[colSums(abs(cor))>0]
    sum(row.names(network.raw)==names(sp.ra2))  #check if matched

    ## 鲁棒性评估robustness simulation
    #input network matrix, percentage of randomly removed species, and ra of all species
    #return the proportion of species remained
    #输入相关矩阵 OTU表格
    Weighted.simu<-rmsimu2(netRaw = network.raw,

                          rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2,
                          abundance.weighted=T,nperm=100)
    head(Weighted.simu)
    Unweighted.simu<-rmsimu2(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2,
                            abundance.weighted=F,nperm=100)
    head(Weighted.simu)
    tem = pst %>% sample_data() %>% .$Group %>% unique() %>% as.character()

    dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),
                     rbind(Weighted.simu,Unweighted.simu),
                     weighted=rep(c("weighted","unweighted"),each=20),
                     Group=tem)
    if (i ==1) {
      dat.f = dat1
    } else if(i != 1){
      dat.f = rbind(dat.f,dat1)
    }

  }

  head(dat.f)

  p = ggplot(dat.f[dat.f$weighted=="weighted",],
             aes(x=Proportion.removed, y=remain.mean, group=Group, color=Group)) +
    geom_line()+
    geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
    xlab("Proportion of species removed")+
    ylab("Proportion of species remained")+
    theme_light()
  p

  p1 = ggplot(dat.f[dat.f$weighted=="unweighted",], aes(x=Proportion.removed,
                                                        y=remain.mean , group=Group, color=Group)) +
    geom_line()+
    geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
    xlab("Proportion of species removed")+
    ylab("Proportion of species remained")+
    theme_light()

  library(patchwork)
  p3 = p|p1

  return(list(p3,dat.f,p,p1))
}












network.efficiency <- function(graph){
  if(is_igraph(graph)==F) warning("Please use a valid iGraph object")
  dd <- 1/shortest.paths(graph)
  diag(dd) <- NA
  efficiency <- mean(dd, na.rm=T)
  #denom <- nrow(dd)*(ncol(dd)-1)
  #sum(dd, na.rm=T)/denom
  return(efficiency)
}

info.centrality.vertex <- function(graph, net=NULL, verbose=F){
  if(is_igraph(graph)==F) warning("Please use a valid iGraph object")
  if(is.null(net)) net <- network.efficiency(graph)
  if(is.numeric(net)==F){
    warning("Please ensure net is a scalar numeric")
    net <- network.efficiency(graph)
  }
  count <- c()
  for(i in 1:length(V(graph))){
    count <- c(count, (net-network.efficiency(igraph::delete.vertices(graph, i)))/net)
    if(verbose){
      print(paste("node",i,"current\ info\ score", count[i], collapse="\t"))
    }
  }
  return(count)
}


info.centrality.network <- function(graph, net=network.efficiency(graph), verbose=F) sum(info.centrality.vertex(graph))
#consider cascade effects: removed species will further influence the remaining nodes
rand.remov2.once<-function(netRaw, rm.num,
                           keystonelist, sp.ra, abundance.weighted=T){
  rm.num2<-ifelse(rm.num > length(keystonelist), length(keystonelist), rm.num)
  id.rm<-sample(keystonelist, rm.num2)
  net.Raw=netRaw %>% as.matrix()
  dim(net.Raw)
  net.new = net.Raw[!names(sp.ra) %in% id.rm, !names(sp.ra) %in% id.rm]   ##remove all the links to these species
  if (nrow(net.new)<2){
    0
  } else {
    sp.ra.new=sp.ra[!names(sp.ra) %in% id.rm]

    if (abundance.weighted){
      net.stength= net.new*sp.ra.new
    } else {
      net.stength= net.new
    }

    sp.meanInteration<-colMeans(net.stength)


    while ( length(sp.meanInteration)>1 & min(sp.meanInteration) <=0){
      id.remain<- which(sp.meanInteration>0)
      net.new=net.new[id.remain,id.remain]
      sp.ra.new=sp.ra.new[id.remain]

      if (abundance.weighted){
        net.stength= net.new*sp.ra.new
      } else {
        net.stength= net.new
      }

      if (length(net.stength)>1){
        sp.meanInteration<-colMeans(net.stength)
      } else{
        sp.meanInteration<-0
      }

    }

    remain.percent<-length(sp.ra.new)/length(sp.ra)

    remain.percent}
}


rmsimuk<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){
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


rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  #-随机挑选出一定百分比的OTU
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw
  #这些节点和其他节点连接全部去除
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    #网络矩阵乘以物种平均丰度，改变相关性值的大小
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  # 每一个节点的平均链接数
  sp.meanInteration<-colMeans(net.stength)

  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.

  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  remain.percent
}

rmsimu2 <-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}


model_compare = function(
    node_table2 = node_table2,
    n = 3,
    padj = FALSE
){

  node_table2$value = 1
  #提取全部的模块及其名称
  module_list = unique(node_table2$group)
  map = node_table2[,2:3] %>% distinct(group, .keep_all = TRUE)
  mytable = node_table2 %>% df_mat(ID,group,value) %>% as.matrix()
  mytable[is.na(mytable)] <- 0

  #--去除少于n个OTU的模块

  if (colSums(mytable)[colSums(mytable)< n] %>% length() == 0) {
    mytable_kp = mytable
  } else {
    mytable_kp = mytable[,-which(colSums(mytable)<n)]
  }


  head(mytable_kp)
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
  i = 1
  for (i in 1:nrow(network_pair)){
    # 全部的需要比对的模块
    module_pair = as.matrix(expand.grid(map_kp$group[which(map_kp$Group==network_pair[i,1])],
                                        map_kp$group[which(map_kp $Group==network_pair[i,2])]))
    overlap = apply(module_pair, 1, FUN= find_overlap, bigtable= mytable_kp)
    only1 = apply(module_pair, 1, FUN= find_only_in_1, bigtable= mytable_kp)
    only2 = apply(module_pair, 1, FUN= find_only_in_2, bigtable= mytable_kp)
    denominator = apply(module_pair, 1, FUN= find_N, mapping=map, bigtable= mytable)
    none = denominator-(overlap + only1 + only2)
    count_table = data.frame(module1 = module_pair[,1],
                             module2 = module_pair[,2], Both=overlap, P1A2=only1, P2A1=only2, A1A2=none)
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

    network1 = network_pair[i,1]
    network2 = network_pair[i,2]
    sig_count = sum(count_table$p_adj<=0.05)
    # count_table$p_adj[3] = 0.001
    if(sig_count>0){
      sig_pairs_table = count_table[which(count_table$p_adj<=0.05),c(1:2)]
      sig_pairs_linked = paste(sig_pairs_table[,1], "-", sig_pairs_table[,2], sep="")
      sig_pairs = paste(sig_pairs_linked, collapse=",")

      sig_pairs_count_table = count_table[which(count_table$p_adj<=0.05),]
      row.names(sig_pairs_count_table) = sig_pairs_linked
      add_one_row = c(network1, network2, sig_count, sig_pairs)
      sig_mod_pairs = rbind(sig_mod_pairs, add_one_row)
      sig_detailed_table = rbind(sig_detailed_table, sig_pairs_count_table)
      sig_detailed_table = sig_detailed_table[-1,]
    }else{
      sig_pairs = "None"
      # add_one_row = c(network1, network2, sig_count, sig_pairs)
      sig_mod_pairs = "none"
      sig_detailed_table = "none"
    }

    print(dim(sig_detailed_table))
    if (i == 1) {
      dat = sig_detailed_table
    } else {
      dat = rbind(dat,sig_detailed_table)
    }
  }


  return(dat)
}


#-utls
find_overlap = function(mods, bigtable){
  vec1 = bigtable[,which(names(bigtable)==mods[1])]
  vec2 = bigtable[,which(names(bigtable)==mods[2])]
  return(sum(vec1*vec2==1))
}

find_only_in_1 = function(mods, bigtable){
  vec1 = bigtable[,which(names(bigtable)==mods[1])]
  vec2 = bigtable[,which(names(bigtable)==mods[2])]
  return(sum(vec1==1 & vec2==0))
}

find_only_in_2 = function(mods, bigtable){
  vec1 = bigtable[,which(names(bigtable)==mods[1])]
  vec2 = bigtable[,which(names(bigtable)==mods[2])]
  return(sum(vec2==1 & vec1==0))
}

find_N = function(mods, mapping, bigtable){
  nwk1 = bigtable[, which(mapping$Group== mapping $Group[which(mapping $group==mods[1])])]
  nwk2 = bigtable[, which(mapping$Group== mapping $Group[which(mapping $group==mods[2])])]
  match_nwk1_nwk2 = sum((rowSums(nwk1) + rowSums(nwk2))>0)
  return(match_nwk1_nwk2)
}

fisher_test = function(x){
  contingency_table <- matrix(unlist(matrix(data.frame(x[3:6]), nrow=2)), nrow=2)
  test_p = fisher.test(contingency_table, alternative = "greater")$p.value
  return(test_p)
}

