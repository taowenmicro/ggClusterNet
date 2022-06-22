#' Microbial related bipartite network analysis
#'
#' @param ps phyloseq Object, contains OTU tables, tax table and map table, represented sequences,phylogenetic tree.
#' @param N filter OTU tables by abundance.The defult, N=0.001, extract the top 0.001 relative abundance of OTU.
#' @param r.threshold The defult, r.threshold=0.6, it represents the correlation that the absolute value
#'  of the correlation threshold is greater than 0.6. the value range of correlation threshold from 0 to 1.
#' @param p.threshold The defult, p.threshold=0.05, it represents significance threshold below 0.05.
#' @param method method for Correlation calculation,method="pearson" is the default value. The alternatives to be passed to cor are "spearman" and "kendall".
#' @param label Whether to add node label.
#' @param group Separate Group.
#' @param env Environmental factor index table which do network analysis with the microbiome.
#' @param envGroup group of env table.
#' @param lay layout which network show.
#' @param path save path of all of network analyse.
#' @param fill fill coulor of node.
#' @param size node size.
#' @param scale Whether relative abundance standardization is required.
#' @param bio Do you need to do a binary network.
#' @param zipi zipi Calculation.
#' @param step Random network sampling times.
#' @param width Save the width of the picture settings.
#' @param height Save the height of the picture setting.
#' @examples
#' path = "./netowrk/"
#' data(ps)
#' ps16s = ps %>% ggClusterNet::scale_micro()
#' psITS = NULL
#' library(phyloseq)
#' ps.merge <- ggClusterNet::merge16S_ITS(ps16s = ps16s,
#'                                        psITS = NULL,
#'                                        N16s = 100)
#' map =  phyloseq::sample_data(ps.merge)
#' head(map)
#' map$Group = "one"
#' phyloseq::sample_data(ps.merge) <- map
#' data(env1)
#' data1 = env1
#' data1$id = row.names(data1)
#' data1 = data1 %>% select(id,everything())
#' envRDA.s = vegan::decostand(env1,"hellinger")
#' data1[,-1] = envRDA.s
#' Gru = data.frame(ID = colnames(data1)[-1],group = "env" )
#' head(Gru)
#' corBionetwork(ps = ps.merge,
#'                             N = 0,
#'                             r.threshold = 0.4, # 相关阈值
#'                             p.threshold = 0.05,
#'                             big = T,
#'                             group = "Group",
#'                             env = data1, # 环境指标表格
#'                             envGroup = Gru,# 环境因子分组文件表格
#'                             # layout = "fruchtermanreingold",
#'                             path = path,# 结果文件存储路径
#'                             fill = "Phylum", # 出图点填充颜色用什么值
#'                             size = "igraph.degree", # 出图点大小用什么数据
#'                             scale = TRUE, # 是否要进行相对丰度标准化
#'                             bio = TRUE, # 是否做二分网络
#'                             zipi = F, # 是否计算ZIPI
#'                             step = 100, # 随机网络抽样的次数
#'                             width = 18,
#'                             label = TRUE,
#'                             height = 10
#' )
#' @return list which contains OTU correlation matrix
#' @author Contact: Tao Wen \email{taowen@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn} yongxin liu \email{yxliu@@genetics.ac.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Tao Wen#, Penghao Xie#, Shengdie Yang, Guoqing Niu, Xiaoyu Liu, Zhexu Ding, Chao Xue, Yong-Xin Liu *, Qirong Shen, Jun Yuan*
#' ggClusterNet: an R package for microbiome network analysis and modularity-based multiple network layouts
#' iMeta 2022,DOI: \url{doi: 10.1002/imt2.32}
#' @export




corBionetwork = function(otu = NULL,
                         tax = NULL,
                         map = NULL,
                         ps = NULL,
                         lab = NULL,
                         N = 0,
                         r.threshold = 0.6, # 相关阈值
                         p.threshold = 0.05,
                         label = FALSE,
                         group = "Group",
                         env = NULL, # 环境指标表格
                         envGroup = NULL,# 环境因子分组文件表格
                         method = "spearman",
                         layout = "fruchtermanreingold",
                         path = "./",# 结果文件存储路径
                         fill = "Phylum", # 出图点填充颜色用什么值
                         size = "igraph.degree", # 出图点大小用什么数据
                         scale = TRUE, # 是否要进行相对丰度标准化
                         bio = TRUE, # 是否做二分网络
                         zipi = FALSE, # 是否计算ZIPI
                         step = 100,
                         width = 20,
                         height = 20,
                         big = TRUE,
                         select_layout = TRUE,
                         layout_net = "model_maptree",
                         clu_method = "cluster_fast_greedy",
                         minsize = 4,
                         maxsize = 14

                         ){
  dir.create(path)
  #--imput data ---------

  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)

  if (scale) {
    ps  = ps %>% ggClusterNet::scale_micro()
  }

  ps_all = ps %>% ggClusterNet::filter_OTU_ps(N)

  mapping = as.data.frame(phyloseq::sample_data(ps))
  mapping$ID = row.names(mapping)
  sample_data(ps) =  mapping

  y = matrix(1,nrow = 14,ncol = length(unique(mapping$Group)))
  layouts = as.character(unique(mapping$Group))

  aa = 1
  plots = list()
  layout = layouts[1]
  layout
  for (layout in layouts) {

    print(layout)
    map <- as.data.frame(phyloseq::sample_data(ps))
    mapsub <- map[map$Group == layout,]
    ps_sub <- ps
    sample_data(ps_sub) <- mapsub
    ps_sub = phyloseq::filter_taxa(ps_sub, function(x) sum(x) > 0 , TRUE)

    if (!is.null(env)) {
      colnames(env)[1] = "ID"
      env_sub <-  env[match(mapsub$ID,env$ID),]
      head(env_sub)
    }

    if (bio) {
      if (!is.null(env)) {

        if (big) {

          result <- corBiostripeBig(data =  env_sub,
                                                  group = envGroup,
                                                  ps = ps_sub,
                                                  r.threshold = r.threshold,
                                                  p.threshold = p.threshold,
                                                  method = method)
        } else {
          result <- corBiostripe(data =  env_sub,
                                               group = envGroup,
                                               ps = ps_sub,
                                               r.threshold = r.threshold,
                                               p.threshold = p.threshold,
                                               method = method)
        }


        #-- extract cor matrix
        occor.r = result[[1]]
        tax = as.data.frame((ggClusterNet::vegan_tax(ps_sub)))

        if (length(tax$filed) != 0) {
          group2 <- data.frame(SampleID = row.names(tax),Group = tax$filed)
        } else {
          group2 <- data.frame(SampleID = row.names(tax),Group = "OTU")
        }

        colnames(envGroup) <-c("SampleID","Group")
        netClu = rbind(envGroup,group2)
        colnames(netClu) <- c("ID","group")
      } else {

        if (big) {
          result <- corBiostripeBig(ps = ps_sub,r.threshold = r.threshold, p.threshold = p.threshold, method = method)
        } else {
          result <- corBiostripe(ps = ps_sub,r.threshold = r.threshold, p.threshold = p.threshold, method = method)
        }

        #-- extract cor matrix
        occor.r = result[[1]]
        tax = as.data.frame((ggClusterNet::vegan_tax(ps_sub)))

        if (length(tax$filed) != 0) {
          group2 <- data.frame(SampleID = row.names(tax),Group = tax$filed)
        } else {
          group2 <- data.frame(SampleID = row.names(tax),Group = "OTU")
        }


        netClu = group2
        colnames(netClu) <- c("ID","group")
      }

    }

    result4 = ggClusterNet::nodeEdge(cor = occor.r)
    igraph  = igraph::graph_from_data_frame(result4[[1]], directed = FALSE, vertices = result4[[2]])

    if (zipi) {
      print("zipi_start")
      #----culculate zi pi
      res = ZiPiPlot(igraph = igraph,method = "cluster_fast_greedy")
      p <- res[[1]]
      ggsave(paste(path,"/",layout,"_ZiPi.pdf",sep = ""),p)
      ZiPi <- res[[2]]
      write.csv(ZiPi ,paste(path,"/",layout,"ZiPi.csv",sep = ""),row.names = FALSE)
    }

    netClu$group = as.factor(netClu$group)


    # result2 = ggClusterNet::PolygonClusterG(cor = occor.r ,nodeGroup = netClu,zoom = 3,zoom2 = 2)

    # library(ggClusterNet)
    if (select_layout) {
      node = NULL
      node = culculate_node_axis(
        cor.matrix = occor.r,
        layout = layout_net,
        seed = 1,
        group = NULL,
        model = FALSE,
        method = clu_method)

    }else if(select_layout) {
      result2 <- model_Gephi.2(cor = cor,
                               method = clu_method,
                               seed = 12
      )
      node = result2[[1]]
    }





    # nodesub = result2[[1]]
    # head(nodesub)

    tax = ps_sub %>%
      ggClusterNet::vegan_tax() %>%
      as.data.frame()

    nodesub1 <- merge(node,tax,by = "row.names",all = T)
    row.names(nodesub1) = nodesub1$Row.names
    nodesub1$Row.names = NULL

    plotcord <- nodesub1 %>%
      dplyr::inner_join(netClu,by =c("elements" = "ID") )


    edges = ggClusterNet::edgeBuild(cor = occor.r,node = node)
    head(edges)
    #-------output---edges and nodes--to Gephi --imput--
    edge_Gephi = data.frame(source = edges$OTU_1,target = edges$OTU_2,correlation =  edges$weight,direct= "undirected",cor =  edges$cor)
    # building node table
    node_Gephi = data.frame(ID= plotcord$elements,plotcord[4:dim(plotcord)[2]],Label = plotcord$elements)
    write.csv(edge_Gephi ,paste(path,"/",layout,"_Gephi_edge.csv",sep = ""),row.names = FALSE)
    write.csv(node_Gephi,paste(path,"/",layout,"_Gephi_node.csv",sep = ""),row.names = FALSE)

    nodepro = ggClusterNet::node_properties(igraph)

    write.csv(nodepro,paste(path,"/",layout,"_node_properties.csv",sep = ""),row.names = FALSE)
    row.names(plotcord) = plotcord$elements
    nodeG = merge(plotcord,nodepro,by = "row.names",all.x = T)
    row.names(nodeG) = nodeG$Row.names
    nodeG$Row.names = NULL
    numna = (dim(nodeG)[2] - 3) : dim(nodeG)[2]
    nodeG[,numna][is.na(nodeG[,numna])] = 0
    head( nodeG)
    p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                    data = edges, size = 0.3,alpha = 0.5) +
      geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = group),pch = 21, data =  nodeG) +
      scale_colour_brewer(palette = "Set1") +
      scale_size(range = c(minsize, maxsize)) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
      labs( title = paste(layout,"network",sep = "_")) + theme_void()
    p0

    head(nodeG)
    tem = nodeG %>% dplyr::filter(elements %in% lab[[1]])

    if (label == TRUE ) {

      if (!is.null(lab)) {
        p1 <- p0 + ggrepel::geom_text_repel(aes(X1, X2,label= elements),size=4, data = tem)
      } else {
        p1 <- p0 + ggrepel::geom_text_repel(aes(X1, X2,label= elements),size=4, data = nodeG)
      }

      plotname = paste(path,"/network_lab_",layout,".pdf",sep = "")
      ggsave(plotname, p1, width = width, height =height)
      # plotname = paste(path,"/network_lab_",layout,".png",sep = "")
      # ggsave(plotname, p1, width = width, height =height)


    }




    p1


    plotname = paste(path,"/network",layout,".pdf",sep = "")
    ggsave(plotname, p0, width = width, height =height)
    # plotname = paste(path,"/network",layout,".png",sep = "")
    # ggsave(plotname, p0, width = width, height =height)
    print("1")
    plots[[aa]] = p0
    rand.g<- erdos.renyi.game(length(V(igraph)), length(E(igraph)),type = c("gnm"))

    data1 = data.frame(network= degree_distribution(igraph, cumulative = FALSE),group = "Erdős–Rényi network",ID = c(1:length(degree_distribution(igraph, cumulative = FALSE))))
    data2 = data.frame(network = degree_distribution(rand.g, cumulative = FALSE) ,group = "network",ID = c(1:length(degree_distribution(rand.g, cumulative = FALSE) )))
    data = rbind(data1,data2)
    p1 <- ggplot(data) +geom_point(aes(x = ID,y = network,group =group,fill = group),pch = 21,size = 2) +
      geom_smooth(aes(x = ID,y = network,group =group,color = group))+
      theme_bw()
    plotname = paste(path,"/Power_law_distribution_",layout,".pdf",sep = "")
    ggsave(plotname, p1, width = width, height =height)

    rand.g.netpro_result<-c()
    for (i in 1:step){
      #####random null model

      rand.g<- erdos.renyi.game(length(V(igraph)), length(E(igraph)),type = c("gnm"))
      tem_netpro_result<- ggClusterNet::net_properties(rand.g)
      rand.g.netpro_result<-cbind(rand.g.netpro_result,tem_netpro_result)
    }
    print("2")
    result_summary<-cbind(rowMeans(rand.g.netpro_result),apply(rand.g.netpro_result,1,sd))
    colnames(result_summary)<-c("Means","SD")


    igraph.weight <- igraph::E(igraph)$weight# 将igraph weight属性赋值到igraph.weight,用于后边做图
    E(igraph)$weight <- NA
    igraph<-igraph::remove.edge.attribute(igraph,"weight")#把边值删除
    netpro_result<- ggClusterNet::net_properties(igraph)
    colnames(netpro_result)<-layout

    sum_net = cbind(netpro_result,result_summary)
    write.csv(sum_net,paste(path,"/",layout,"_net_VS_erdos_properties.csv",sep = ""),row.names = TRUE)
    print("3")
    y = as.data.frame(y)
    colnames(y) = layouts
    # head(y)
    y[layout] = netpro_result[,1]
    row.names(y) = row.names(netpro_result)
    aa = aa+1
  }

  plotname = paste(path,"/network_all.pdf",sep = "")
  p  = ggpubr::ggarrange(plotlist = plots, common.legend = TRUE, legend="right")
  return(list(p,y,edges,nodeG,occor.r))

}



#' Construct a network layout. Calculate the layout according to grouping and random distribution
#'
#' @param cor Correlation matrix
#' @param nodeGroup Classification information of network nodes
#' @param zoom Set the distance between modules
#' @param zoom2 Scaling module radius size
#' @examples
#' data
#' data(ps)
#' result = corMicro (ps = ps,N = 0.02,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # building the node group
#' ps_net = result[[3]]
#' vegan_tax <-  function(physeq){
#' tax <-  tax_table(physeq)
#'
#' return(as(tax,"matrix"))
#' }
#' tax_table = as.data.frame(vegan_tax(ps_net))
#' group = as.data.frame(tax_table)
#' group$ID = row.names(group)
#' netClu = data.frame(ID = row.names(group),group = group$Phylum)
#' netClu$group = as.factor(netClu$group)
#' result2 = PolygonRrClusterG (cor = cor,nodeGroup =netClu )
#' node = result2[[1]]
#'
#'
#' @return result2 Which contains 2 lists.Result2[[1]], consists of OTU and its corresponding coordinates.
#' result2[[2]], consists of the network center coordinates of each group
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export



PolygonRrClusterG = function(cor = cor,nodeGroup =netClu,zoom = 1,zoom2 = 1,bio = F){
  num = length(levels(nodeGroup$group))
  xs = as.data.frame(table(nodeGroup$group))
  r = xs$Freq/10 *zoom

  # Calculate angle according to group
  arg = seq(0,360,360/(length(r))) - 180
  # i = 1
  rsum = sum(r)
  x= rep(0,length(r))
  y = rep(0,length(r))
  for (i in 1:length(r)) {
    x[i] = (rsum + r[i])* sin(arg[i]* 3.14/180)
    y[i] = (rsum + r[i])* cos(arg[i]* 3.14/180)
  }
  if (bio == T) {
    for (i in 1:length(r)) {
      x[i] = 0
      y[i] = 0
    }
  }
  da = data.frame(x = x,y = y)
  for (i in 1:length(levels(nodeGroup$group))) {
    # Extract all otu in this group
    as = dplyr::filter(nodeGroup, group == levels(nodeGroup$group)[i])
    if (length(as$ID) == 1) {
      data = cbind(da[i,1],da[i,2] )
      data =as.data.frame(data)
      row.names(data ) = as$ID
      data$elements = row.names(data )
      colnames(data)[1:2] = c("X1","X2")
    }
    as$ID = as.character(as$ID)

    if (length(as$ID)!=1 ) {
      m = cor[as$ID,as$ID]

      d  =m
      d <- sna::as.edgelist.sna(d)
      # if (is.list(d))
      # d <- d[[1]]
      n <- attr(d, "n")
      # 提取半径
      s = r[i]
      s = s * zoom2

      data = cbind(sin(2 * pi * ((0:(n - 1))/n))*s +da[i,1], cos(2 * pi * ((0:(n - 1))/n))*s +da[i,2])

      data =as.data.frame(data)
      row.names(data ) = row.names(m)
      data$elements = row.names(data )
      colnames(data)[1:2] = c("X1","X2")

    }

    if (i == 1) {
      oridata = data
    }
    if (i != 1) {
      oridata = rbind(oridata,data)
    }

  }
  plotcord = oridata[match(oridata$elements,row.names(cor )),]

  return(list(plotcord,da))
}

