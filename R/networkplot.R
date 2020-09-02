#' Microbial related network
#'
#' @param ps phyloseq Object, contains OTU tables, tax table and map table, represented sequences,phylogenetic tree.
#' @param N filter OTU tables by abundance.The defult, N=0.001, extract the top 0.001 relative abundance of OTU.
#' @param r.threshold The defult, r.threshold=0.6, it represents the correlation that the absolute value
#'  of the correlation threshold is greater than 0.6. the value range of correlation threshold from 0 to 1.
#' @param p.threshold The defult, p.threshold=0.05, it represents significance threshold below 0.05.
#' @param method method for Correlation calculation,method="pearson" is the default value. The alternatives to be passed to cor are "spearman" and "kendall".
#' @param label Whether to add node label.
#' @param group Separate Group.
#' @param lay layout which network show
#' @param path save path of all of network analyse.
#' @param fill fill coulor of node
#' @param size node size
#' @param zipi zipi Calculation
#' @param step Random network sampling times
#' @examples
#' data(ps)
#' result = network (ps = ps,N = 0.001,r.threshold=0.6,p.threshold=0.05,label = FALSE,path = path ,zipi = TRUE)
#' result[[1]]
#' result[[2]]
#' @return list which contains OTU correlation matrix
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export



network = function(otu = NULL,
                   tax = NULL,
                   map = NULL,
                   ps = NULL,
                   N = 0.001,
                   r.threshold = 0.6,
                   p.threshold = 0.05,
                   method = "pearson",
                   label = FALSE,
                   group = "Group",
                   lay = "fruchtermanreingold",
                   path = "./",
                   fill = "Phylum",
                   size = "igraph.degree",
                   scale = TRUE,
                   zipi = FALSE,
                   clu_method = "cluster_fast_greedy",
                   step = 100){


  #--imput data ---------
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  #---------------------------data washing-------------------------------------------------
  #transform to relative abundance
  if (scale == TRUE) {
    ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )

  } else {
    ps_rela <- ps
  }

  ps_sub = phyloseq::filter_taxa(ps_rela, function(x) mean(x) > N , TRUE)#select OTUs according to  relative abundance


  # extract map table
  design = mapping= as.data.frame(phyloseq::sample_data(ps_sub))
  colnames(mapping) <- "Group"
  otu_table = as.data.frame(t(vegan_otu(ps_sub)))
  tax_table = as.data.frame(vegan_tax(ps_sub ))

  # network = otu_table
  y = matrix(1,nrow = 14,ncol = length(unique(mapping$Group)))
  #--transmit N
  d = N


  layouts = as.character(unique(mapping$Group))
  mapping$ID = row.names(mapping)
  ##################---------------------------------------calculate network---------------------------------------------------
  aa = 1
  plots = list()
  # layout = layouts[1]

  for (layout in layouts) {

    # xx = dplyr::filter(as.tibble(mapping), Group %in% layout)
    xx <- mapping[mapping$Group ==  layout,]

    network_sub =  otu_table[xx$ID]
    n=ncol(network_sub)
    network_sub[n+1]=apply(network_sub[c(1:nrow(network_sub)),],1,sum)
    network_sub=network_sub[network_sub[n+1] > 0,1:n]

    print(layout)
    occor = psych::corr.test(t(network_sub),use="pairwise",method= method,adjust="fdr",alpha=.05)
    print("over")
    occor.r = occor$r #
    occor.p = occor$p #
    occor.r[occor.p > p.threshold|abs(occor.r)< r.threshold] = 0

    # if (bio == TRUE) {
    #   # intersect(row.names(occor.r),as.character(row.names(tax_table[as.character(tax_table$tax3)== "ARG",])))
    #   a <- row.names(occor.r) %in% intersect(row.names(occor.r),as.character(row.names(tax_table[as.character(tax_table$tax3)== "ARG",])))
    #   a
    #   occor.r[a,a] = 0
    #   a <- row.names(occor.r) %in% intersect(row.names(occor.r),as.character(row.names(tax_table[as.character(tax_table$tax3)== "MEG",])))
    #   a
    #   occor.r[a,a] = 0
    #
    # }
    result4 = nodeEdge(cor = occor.r)
    # extract edge file
    edge = result4[[1]]
    # extract node file
    node = result4[[2]]
    #--
    igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
    # nodepro = node_properties(igraph)
    if (zipi == TRUE) {
      #----culculate zi pi
      res = ZiPiPlot(igraph = igraph,method = clu_method)
      p <- res[[1]]
      ggsave(paste(path,"/",layout,"_ZiPi.pdf",sep = ""),p)
      ZiPi <- res[[2]]
      write.csv(ZiPi ,paste(path,"/",layout,"ZiPi.csv",sep = ""),row.names = FALSE)
    }

    g <- network::network(occor.r, directed=FALSE)
    g  = g
    m <- network::as.matrix.network.adjacency(g)  # get sociomatrix

    # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
    # plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))

    plotcord <- selectlayout(m = m,layout = lay)


    colnames(plotcord) = c("X1", "X2")
    plotcord$elements <- colnames(occor.r)
    edglist <- network::as.matrix.network.edgelist(g)
    edglist = as.data.frame(edglist)

    network::set.edge.value(g,"weigt",occor.r)
    edglist$weight = as.numeric(network::get.edge.attribute(g,"weigt"))
    edges <- data.frame(plotcord[edglist[, 1], ], plotcord[edglist[, 2], ])
    edges$weight = as.numeric(network::get.edge.attribute(g,"weigt"))
    aaa = rep("a",length(edges$weight))
    for (i in 1:length(edges$weight)) {
      if (edges$weight[i]> 0) {
        aaa[i] = "+"
      }
      if (edges$weight[i]< 0) {
        aaa[i] = "-"
      }
    }
    #add to edge table
    edges$wei_label = as.factor(aaa)
    colnames(edges) <- c("X1", "Y1","OTU_1", "X2", "Y2","OTU_2","weight","wei_label")
    edges$midX <- (edges$X1 + edges$X2)/2
    edges$midY <- (edges$Y1 + edges$Y2)/2

    ##plotcord add tax
    row.names(plotcord) = plotcord$elements
    network_tax =tax_table
    res = merge(plotcord,network_tax,by = "row.names",all = F)

    row.names(res) = res$Row.names
    res$Row.names = NULL
    plotcord = res

    xx = data.frame(mean  =rowMeans(otu_table))
    plotcord = merge(plotcord,xx,by = "row.names",all = FALSE)
    row.names(plotcord) = plotcord$Row.names
    plotcord$Row.names = NULL

    #-------output---edges and nodes--to Gephi --imput--
    edge_Gephi = data.frame(source = edges$OTU_1,target = edges$OTU_2,correlation =  edges$weight,direct= "undirected",cor =  edges$wei_label)
    # building node table
    node_Gephi = data.frame(ID= plotcord$elements,plotcord[4:dim(plotcord)[2]],Label = plotcord$elements)

    write.csv(edge_Gephi ,paste(path,"/",layout,"_Gephi_edge.csv",sep = ""),row.names = FALSE)
    write.csv(node_Gephi,paste(path,"/",layout,"_Gephi_node.csv",sep = ""),row.names = FALSE)

    nodepro = node_properties(igraph)

    write.csv(nodepro,paste(path,"/",layout,"_node_properties.csv",sep = ""),row.names = FALSE)

    nodeG = merge(plotcord,nodepro,by = "row.names",all  =FALSE)
    head(nodeG)
    row.names(nodeG) = nodeG$Row.names
    nodeG$Row.names = NULL
    #--
    plotnode <- nodeG %>% dplyr::select(c("X1" , "X2","elements",fill, "mean",size))
    colnames(plotnode) <- gsub(fill,"XXXX",colnames(plotnode))
    colnames(plotnode) <- gsub(size,"YYYY",colnames(plotnode))

    pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(wei_label)),
                                    data = edges, size = 0.5) +
      geom_point(aes(x = X1, y = X2,size = YYYY,fill =XXXX),pch = 21, data =  plotnode) + scale_colour_brewer(palette = "Set1") +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
      labs( title = paste(layout,"network",sep = "_")) + theme_void()

    if (label == TRUE ) {
      pnet <- pnet +  geom_text_repel(aes(X1, X2,label=XXXX),size=4, data = plotnode)
    }

    plotname = paste(path,"/network",layout,".pdf",sep = "")
    ggsave(plotname, pnet, width = 12, height =8)

    plots[[aa]] = pnet
    #----
    rand.g<- erdos.renyi.game(length(V(igraph)), length(E(igraph)),type = c("gnm"))
    # degree_distribution

    ## We compare the degree of random network and this network
    data1 = data.frame(network= degree_distribution(igraph, cumulative = FALSE),group = "E–R network",ID = c(1:length(degree_distribution(igraph, cumulative = FALSE))))
    data2 = data.frame(network = degree_distribution(rand.g, cumulative = FALSE) ,group = "network",ID = c(1:length(degree_distribution(rand.g, cumulative = FALSE) )))
    data = rbind(data1,data2)
    p1 <- ggplot(data) +geom_point(aes(x = ID,y = network,group =group,fill = group),pch = 21,size = 2) +
      geom_smooth(aes(x = ID,y = network,group =group,color = group))+
      theme_bw()
    plotname = paste(path,"/Power_law_distribution_",layout,".pdf",sep = "")
    ggsave(plotname, p1, width = 8, height =6)

    #---------------------

    rand.g.netpro_result<-c()
    for (i in 1:step){
      #####random null model
      rand.g<- erdos.renyi.game(length(V(igraph)), length(E(igraph)),type = c("gnm"))
      tem_netpro_result<-net_properties(rand.g)
      rand.g.netpro_result<-cbind(rand.g.netpro_result,tem_netpro_result)
    }
    #--------计算step次的结果
    # rand.g.netpro_result

    # --简化结果--------对随机矩阵结果求均值和sd值
    result_summary<-cbind(rowMeans(rand.g.netpro_result),apply(rand.g.netpro_result,1,sd))
    colnames(result_summary)<-c("Means","SD")
    head(result_summary)


    ###网络边的赋值及其设置
    igraph.weight <- E(igraph)$weight# 将igraph weight属性赋值到igraph.weight,用于后边做图
    E(igraph)$weight <- NA
    igraph<-remove.edge.attribute(igraph,"weight")#把边值删除
    netpro_result<- net_properties(igraph)
    colnames(netpro_result)<-layout

    sum_net = cbind(netpro_result,result_summary)
    write.csv(sum_net,paste(path,"/",layout,"_net_VS_erdos_properties.csv",sep = ""),row.names = TRUE)

    y = as.data.frame(y)
    colnames(y) = layouts
    # head(y)
    y[layout] = netpro_result[,1]
    row.names(y) = row.names(netpro_result)
    aa = aa+1
  }
  plotname = paste(path,"/network_all.pdf",sep = "")
  p  = ggpubr::ggarrange(plotlist = plots, common.legend = TRUE, legend="right")
  return(list(p,y))
}




selectlayout <- function(m,layout = "fruchtermanreingold"){
  if (layout == "fruchtermanreingold") {
    plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))

  }

  if (layout == "circle") {
    plotcord <-  data.frame(sna::gplot.layout.circle(m, NULL))
  }

  if (layout == "kamadakawai") {
    plotcord <-  data.frame(sna::gplot.layout.adj(m, NULL))
  }
  if (layout == "adj") {
    plotcord <-  data.frame(sna::gplot.layout.kamadakawai(m, NULL))
  }
  if (layout == "circrand") {
    plotcord <-  data.frame(sna::gplot.layout.circrand(m, NULL))
  }
  if (layout == "eigen") {
    plotcord <-  data.frame(sna::gplot.layout.eigen(m, NULL))
  }
  if (layout == "geodist") {
    plotcord <-  data.frame(sna::gplot.layout.geodist(m, NULL))
  }
  if (layout == "hall") {
    plotcord <-  data.frame(sna::gplot.layout.hall(m, NULL))
  }

  if (layout == "mds") {
    plotcord <-  data.frame(sna::gplot.layout.mds(m, NULL))
  }

  if (layout == "princoord") {
    plotcord <-  data.frame(sna::gplot.layout.princoord(m, NULL))
  }
  if (layout == "random") {
    plotcord <-  data.frame(sna::gplot.layout.random(m, NULL))
  }
  if (layout == "rmds") {
    plotcord <-  data.frame(sna::gplot.layout.rmds(m, NULL))
  }
  if (layout == "segeo") {
    plotcord <-  data.frame(sna::gplot.layout.segeo(m, NULL))
  }
  if (layout == "seham") {
    plotcord <-  data.frame(sna::gplot.layout.seham(m, NULL))
  }
  if (layout == "spring") {
    plotcord <-  data.frame(sna::gplot.layout.spring(m, NULL))
  }
  if (layout == "springrepulse") {
    plotcord <-  data.frame(sna::gplot.layout.springrepulse(m, NULL))
  }
  if (layout == "target") {
    plotcord <-  data.frame(sna::gplot.layout.target(m, NULL))
  }

  return(plotcord)
}




