natural.con.microp = function(
    ps = ps,
    Top = 500,
    corg = NULL,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman",
    norm = TRUE,
    end = 400,
    start = 0,
    con.method = "pulsar"

){

  # otutab<- ps %>%
  #   vegan_otu() %>%
  #   as.data.frame()
  # dim(otutab)

  # id <- sample_data(ps)$Group %>% unique()
  if (!is.null(corg)) {
    id = names(corg)
  } else if (is.null(corg)){

  }

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



    dat = natural.con.micro(cor = cor,start = start,
                            end = end,
                            norm = norm,
                            method = con.method
    )
    dat$Group = as.character(id[i])
    if (i ==1) {
      datall = dat
    }else if (i != 1){
      datall = rbind(datall,dat)
    }

  }


  datall$Natural.connectivity = Re(datall$Natural.connectivity)
  p=ggplot(datall, aes(Num.of.remove.nodes, Natural.connectivity,color=Group)) +
    geom_point(alpha = 0.1) +
    geom_smooth(level = F) +
    theme_set(theme_bw(base_size = 25))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black")) + theme_bw()
  # p
  return(list(p,datall))
}


natural.con.micro  =function(cor,
                             start,
                             end,
                             norm = F,
                             method = "pulsar"# "xph"
){
  A = c()
  B = c()
  j = 1
  for (j in start:end) {
    num = length(row.names(cor)) - j
    tem = sample(row.names(cor),num)
    cor.z = cor[tem,tem] %>% as.matrix()
    # cor.z[1:50,1:5]
    # igraph = make_igraph(cor)
    #存在某些情况计算不出来相关系数，定义相关为0
    cor[is.na(cor.z)]<-0
    cor.z[is.na(cor.z)] = 0
    #-去除自相关的点
    diag(cor.z) <- 0
    #-查看网络边的数量
    # tem = abs(cor.z)>0 %>% as.vector()
    # tem2 = as.matrix(tem) %>% as.vector() %>% table() %>% as.data.frame()
    # tem2[2,2]/2


    sum(tem>0)/2
    #网络中节点的数量
    sum(colSums(abs(cor.z))>0)  # node number: number of species with at least one linkage with others.
    # #去除没有任何相关的节点.
    # network.raw<-cor[colSums(abs(cor.z))>0,colSums(abs(cor.z))>0]
    A[j] = j
    if (method == "pulsar" ) {
      B[j] = natural.connectivity(G = cor.z,eig =NULL,norm = norm)
    } else if(method == "xph"){
      B[j] =nc(cor.z)

    }
  }







  dat = data.frame(Num.of.remove.nodes = A,Natural.connectivity = B)
  return(dat)
}




nc <- function(adj_matrix = cor) {
  #获取 0-1 矩阵，1 表示节点间存在边，0 表示不存在边
  adj_matrix <- as.matrix(adj_matrix)
  adj_matrix[abs(adj_matrix) != 0] <- 1
  diag(adj_matrix)<-0
  #矩阵的特征分解，获取特征值 λ
  lambda <- eigen(adj_matrix, only.values = TRUE)$values
  lambda <- sort(lambda, decreasing = TRUE)

  #计算“平均特征根”，获得自然连通度
  lambda_sum <- 0
  N = length(lambda)
  for (i in 1:N) lambda_sum = lambda_sum + exp(lambda[i])
  lambda_average <- log(lambda_sum/N, base = exp(1))
  lambda_average
}

