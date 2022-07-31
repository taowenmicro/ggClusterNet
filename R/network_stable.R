
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

rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
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

