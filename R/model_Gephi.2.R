#' # The algorithm imitate Gephi's network layout. while, it saves time and the calculation speed is faster,This is the upgraded version.
#'
#' @title This is the upgraded version.model_Gephi.2:The algorithm imitate Gephi's network layout. while, it saves time and the calculation speed is faster
#' @description Enter correlation matrix, calculate network modules, and Calculate the coordinates of the node.
#' @param cor Correlation matrix Clustering Algorithm
#' @param method Clustering Algorithm
#' @param seed Set random seed
#' @details
#' By default, returns a list
#' \itemize{
#' \item{cluster_fast_greedy: }
#' \item{cluster_walktrap: }
#' \item{cluster_edge_betweenness: }
#' \item{cluster_spinglass: }
#' }
#' @examples
#' data(ps)
#' result = corMicro (ps = ps,N = 100,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # building the node group
#' result2 <- model_Gephi.2(cor = cor,
#' method = "cluster_fast_greedy",
#' seed = 12
#' )
#' node = result2[[1]]
#' node = result2[[2]]
#' @return list
#' @author Contact: Tao Wen \email{taowen@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn} yongxin liu \email{yxliu@@genetics.ac.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Tao Wen#, Penghao Xie#, Shengdie Yang, Guoqing Niu, Xiaoyu Liu, Zhexu Ding, Chao Xue, Yong-Xin Liu *, Qirong Shen, Jun Yuan*
#' ggClusterNet: an R package for microbiome network analysis and modularity-based multiple network layouts
#' iMeta 2022,DOI: \url{doi: 10.1002/imt2.32}
#' @export



model_Gephi.2 <- function(cor = cor,method = "cluster_fast_greedy",seed = 2){
  #---节点的模块化计算
  corr = cor
  # Construct Edge File
  edges <- data.frame(from = rep(row.names(corr), ncol(corr)),
                      to = rep(colnames(corr), each = nrow(corr)),
                      r = as.vector(corr)
  )
  # Extract half of the matrix, and remove the diagonal correlation (self related to yourself)
  edges <- dplyr::filter(edges, as.vector(lower.tri(corr)))
  colnames(edges)[3] = "weight"
  #---Set the sign of the edge
  # E.color <- edges$weight
  edges$direction <- ifelse(edges$weight>0, "pp",ifelse(edges$weight<0, "np","ns"))
  node = data.frame(name = unique(c(as.character(edges$from),as.character( edges$to))))
  row.names(node) = node$name

  # Output igraph object
  igraph  = igraph::graph_from_data_frame(edges, directed = FALSE, vertices = node)
  if (method == "cluster_walktrap" ) {
    fc <- cluster_walktrap(igraph,weights =  abs(E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }

  if (method == "cluster_edge_betweenness" ) {
    fc <- cluster_edge_betweenness(igraph,weights =  abs(E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }
  if (method == "cluster_fast_greedy" ) {
    fc <- cluster_fast_greedy(igraph,weights =  abs(E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }
  if (method == "cluster_spinglass" ) {
    fc <- cluster_spinglass(igraph,weights =  abs(E(igraph)$weight))# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }

  modularity <- modularity(igraph,membership(fc))
  # Modularity
  modularity
  #-Extraction module

  netClu = data.frame(ID = names(membership(fc)),group = as.vector(membership(fc)))
  dim(netClu)

  table(netClu$group)

  result4 = nodeEdge(cor = cor)
  edge = result4[[1]]
  node = result4[[2]]
  #--
  igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
  igraph.degree<-igraph::degree(igraph) %>% as.data.frame()
  colnames(igraph.degree) = "degree"
  igraph.degree$ID = row.names(igraph.degree)
  dim(igraph.degree)
  netClu <- netClu %>%
    dplyr::left_join(igraph.degree,na_matches = "never")
  netClu$degree[is.na(netClu$degree)] = 0
  netClu <- netClu %>%
    dplyr::arrange(desc(degree))
  sumtav <- netClu %>% dplyr::group_by(group) %>%
    dplyr::summarise(sum(degree))
  colnames(sumtav) = c("group","degree")
  tab0 <- sumtav %>%
    dplyr::arrange(desc(degree))
  tab0$group = as.character(tab0$group)

  tab1 = as.data.frame(table(netClu$group)) %>%
    dplyr::arrange(desc(Freq))
  colnames(tab1)[1] = "group"
  tab1$group = as.character(tab1$group)
  tab3 <- tab0 %>%
    dplyr::left_join(tab1,by = "group")

  num.node <- dim(cor)[1]

  for (N in 1: num.node) {
    A = 1 + (7*(N + 1)*N )/2 - N
    if (A >= num.node) {

     break
    }
    n = N - 1
    print(n)
  }


  # n = (sqrt((num.node-1)/3) - 1) %>% floor()
  wai.mode = num.node - (1 + (7*(n + 1)*n )/2 - n)
  dat = data.frame(x = 0,y = 0)
  i = 1
  for (i in 1:n) {
    t <- seq(0, 2*pi, length.out = 7*i)
    t = t[-1]
    x <- sin(t)*i
    y <- cos(t)*i
    add = data.frame(x = x,y = y)
    dat = rbind(dat,add)

    if (i== n) {
      i = i + 1
      t <- seq(0, 2*pi, length.out = (wai.mode + 1))
      t = t[-1]
      x <- sin(t)*i
      y <- cos(t)*i
      add = data.frame(x = x,y = y)
      dat = rbind(dat,add)
    }

  }
  row.names(dat) = row.names(cor)
  dat$elements = row.names(cor)
  colnames(dat)[1:2] = c("X1","X2")

  dat0 = dat

  # head(dat)
  # j  = 1
  axis.node = c()

  for (j in 1:dim(tab3)[1]) {

    if (dim(dat)[1] <= 2 | tab3$Freq[j] == tab3$Freq[dim(tab3)[1]]) {
      lacat = row.names(dat) [1:tab3$Freq[j]]
    }

    if (dim(dat)[1] > 2 & tab3$Freq[j] != tab3$Freq[dim(tab3)[1]]) {
      set.seed(seed)
      axis_mod = sample(1:dim(dat)[1],100,replace = TRUE) %>% sort()

      d <- dist(dat[,-3]) %>% as.matrix()
      id = dat[axis_mod[1],]$elements
      id2 = d[id,] %>% sort()
      lacat = c(names(id2[1:(tab3$Freq[j])]) )
    }

    new.dat = dat[lacat,]
    New.locat  = netClu$ID[netClu$group == as.numeric(tab3$group[j])]
    row.names(new.dat) = New.locat
    new.dat$elements = New.locat

    # ggplot(new.dat,aes(x = X1,y = X2)) + geom_point() + geom_text(aes(label= elements))

    if (j == 1) {
      axis.node = new.dat
    }
    if (j != 1) {
      axis.node = rbind(axis.node,new.dat)
    }
    dat = dat[dat$elements %in% lacat == FALSE,]

  }

  return(list(axis.node,dat0,netClu))


}

