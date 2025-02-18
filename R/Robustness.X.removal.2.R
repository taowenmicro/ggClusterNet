
Robustness.Random.removal.2 = function(
    ps = ps,
    corg = NULL,
    Top = 500,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman"
){
  otutab<- ps %>%
    vegan_otu() %>%
    as.data.frame()
  dim(otutab)
  
  if (!is.null(corg)) {
    id = names(corg)
  } else if (is.null(corg)){
    
  }
  
  id <- sample_data(ps)$Group %>% unique()
  i  = 1
  #计算每个物种的平均丰度，使用测序深度标准化
  sp.ra<-colMeans(otutab)/mean(rowSums(otutab))   #relative abundance of each species
  
  # library(ggClusterNet)
  for (i in 1:length(id)){
    
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
                           abundance.weighted=TRUE,nperm=100)
    head(Weighted.simu)
    Unweighted.simu<-rmsimu2(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2,
                             abundance.weighted=FALSE,nperm=100)
    head(Weighted.simu)
    
    tem = ps %>%
      scale_micro() %>%
      subset_samples.wt("Group", c(id[i])) %>%
      subset_taxa.wt("OTU",row.names(cor)) %>%
      # filter_OTU_ps(Top) %>%
      sample_data() %>%
      .$Group %>%
      unique() %>% as.character()
    
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









Robustness.Targeted.removal.2 = function(
    ps = ps,
    corg = NULL,
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
  
  if (!is.null(corg)) {
    id = names(corg)
    id <- phyloseq::sample_data(ps)$Group %>% unique()
  } else if (is.null(corg)){
    id <- phyloseq::sample_data(ps)$Group %>% unique()
    
  }
  
  
  #计算每个物种的平均丰度，使用测序深度标准化
  sp.ra<-colMeans(otutab)/mean(rowSums(otutab))   #relative abundance of each species
  
  for (i in 1:length(id)){
    
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
    
    otutab<- ps %>%
      scale_micro() %>%
      subset_samples.wt("Group", c(id[i])) %>%
      # filter_OTU_ps(Top) %>%
      subset_taxa.wt("OTU",row.names(cor)) %>%
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
      model = res[[2]] %>%  dplyr::filter(roles == "Module hubs")
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
        dplyr::filter(!is.na(igraph.degree) ) %>%
        dplyr::arrange(desc(igraph.degree))
      head(ret3)
      tem = round(length(ret3$igraph.degree) * 0.05,0)
      module.hub = row.names(ret3)[1:tem]
    }
    
    
    rmsimu<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=TRUE,nperm=100){
      
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
                          abundance.weighted=TRUE,nperm=100)
    Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub),
                            keystonelist=module.hub,
                            sp.ra=sp.ra2, abundance.weighted=FALSE,nperm=100)
    
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








