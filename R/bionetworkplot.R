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

corBionetwork = function(otu = NULL,
                   tax = NULL,
                   map = NULL,
                   ps = NULL,
                   N = 0.001,
                   r.threshold = 0.6,
                   p.threshold = 0.05,
                   label = FALSE,
                   group = "Group",
                   env = NULL,
                   envGroup = NULL,
                   method = "spearman",
                   layout = "fruchtermanreingold",
                   path = "./",
                   fill = "Phylum",
                   size = "igraph.degree",
                   scale = TRUE,
                   bio = FALSE,
                   zipi = FALSE,
                   step = 100,
                   width = 8,
                   height = 6,
                   allnode = TRUE,
                   ncol = 3,
                   nrow = 1){
  #--imput data ---------
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  #---------------------------data washing-------------------------------------------------
  #transform to relative abundance
  if (scale == TRUE) {
    ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
  } else {
    ps_rela <- ps
  }
  ps_all = phyloseq::filter_taxa(ps_rela, function(x) mean(x) > N , TRUE)#select OTUs according to  relative abundance
  # extract map table otu tax table-------------------
  mapping = as.data.frame(phyloseq::sample_data(ps_all))
  mapping$ID = row.names(mapping)
  # colnames(mapping)[1] <- "ID"
  sample_data(ps_all) =  mapping
  y = matrix(1,nrow = 14,ncol = length(unique(mapping$Group)))
  layouts = as.character(unique(mapping$Group))

  ##################---------------------------------------calculate network---------------------------------------------------
  aa = 1
  plots = list()
  # layout = layouts[1]
  for (layout in layouts) {

    print(layout)
    map <- as.data.frame(phyloseq::sample_data(ps_all))
    map$Group = as.character(map$Group)
    mapsub <- map[map$Group == layout,]
    ps_sub <- ps_all
    sample_data(ps_sub) <- mapsub
    ps_sub = filter_taxa(ps_sub, function(x) sum(x) > 0 , TRUE)# remove 0 of all sample


    # match(as.character(mapsub$ID),colnames(env))

    # colnames(env)[1] = "ID"
    # xx = dplyr::filter(as.tibble(map), Group %in% layout)

    if (length(match(row.names(env),mapsub$ID)[!is.na(match(row.names(env),mapsub$ID))]) == 0) {
     env = t(env)
     env = as.data.frame(env)
    }

    env$ID = row.names(env)

    env_sub <-  env[match(as.character(mapsub$ID),as.character(env$ID)),]
    # env_sub <- env_sub[,c(ncol(env_sub),1:c(ncol(env_sub) - 1))]
    env_sub <- env_sub %>% select(ID,everything())

    if (bio == TRUE) {

      result <- corBiostripe(data =  env_sub,group = envGroup,ps = ps_sub,r.threshold = r.threshold, p.threshold = p.threshold, method = method)
      #-- extract cor matrix
      occor.r = result[[1]]
      tax = as.data.frame((vegan_tax(ps_sub)))
      if (length(tax$filed) != 0) {
        group2 <- data.frame(SampleID = row.names(tax),Group = tax$filed)
      } else {
        group2 <- data.frame(SampleID = row.names(tax),Group = "OTU")
      }
      colnames(envGroup) <-c("SampleID","Group")
      row.names(envGroup) =  envGroup$SampleID
      netClu = rbind(group2,envGroup)
      colnames(netClu) <- c("ID","group")
    }

    result4 = nodeEdge(cor = occor.r)

    # extract edge file
    # extract node file
    #--
    igraph  = igraph::graph_from_data_frame(result4[[1]], directed = FALSE, vertices = result4[[2]])
    # nodepro = node_properties(igraph)
    print("igraph_over")
    if (zipi == TRUE) {
      print("zipi_start")
      #----culculate zi pi
      res = ZiPiPlot(igraph = igraph,method = "cluster_fast_greedy")
      p <- res[[1]]
      ggsave(paste(path,"/",layout,"_ZiPi.pdf",sep = ""),p)
      ZiPi <- res[[2]]
      write.csv(ZiPi ,paste(path,"/",layout,"ZiPi.csv",sep = ""),row.names = FALSE)
    }

    netClu$group = as.factor(netClu$group)

    result2 = PolygonClusterG (cor = occor.r ,nodeGroup = netClu,zoom = 3,zoom2 = 2 )


    nodesub = result2[[1]]
    print("PolygonRrClusterG")

    ### nodeadd 节点注释的简单封装，便捷实用otu表格和分组文件进行注释
    plotcord <- nodesub %>%
      inner_join(netClu,by =c("elements" = "ID") )
    #-----计算边#--------
    print("start_edgebuilding")
    edges = edgeBuild(cor = occor.r,plotcord = nodesub)
    edges = edges %>% filter(wei_label != "a")
    head(plotcord)
    print(dim(edges))
    #-------output---edges and nodes--to Gephi --imput--
    edge_Gephi = data.frame(source = edges$OTU_1,target = edges$OTU_2,correlation =  edges$weight,direct= "undirected",cor =  edges$wei_label)
    # building node table
    node_Gephi = data.frame(ID= plotcord$elements,plotcord[4:dim(plotcord)[2]],Label = plotcord$elements)

    write.csv(edge_Gephi ,paste(path,"/",layout,"_Gephi_edge.csv",sep = ""),row.names = FALSE)
    write.csv(node_Gephi,paste(path,"/",layout,"_Gephi_node.csv",sep = ""),row.names = FALSE)
    print("gephi_over")

    nodepro = node_properties(igraph)

    write.csv(nodepro,paste(path,"/",layout,"_node_properties.csv",sep = ""),row.names = FALSE)
    row.names(plotcord) = plotcord$elements

    if (allnode == TRUE) {
      nodeG = merge(plotcord,nodepro,by = "row.names",all.x = TRUE)
    }else {
      nodeG = merge(plotcord,nodepro,by = "row.names",all = FALSE)
    }
    row.names(nodeG) = nodeG$Row.names
    nodeG$Row.names = NULL
    nodeG[is.na(nodeG)] = 0

    # if (fill == group) {
    #   fill = "group"
    # }

    head(nodeG)
    plotnode <- nodeG
    # plotnode <- nodeG %>% dplyr::select(c("X1" , "X2","elements",fill, size))
    # colnames(plotnode) <- gsub(fill,"XXXX",colnames(plotnode))
    # colnames(plotnode) <- gsub(size,"YYYY",colnames(plotnode))
    # head(plotnode)
    # !!sym(a)
    # fill = "group"
    # size = "igraph.betweenness"
    lab = "elements"
    ### 出图



    pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(wei_label)),
                                    data = edges, size = 0.5) +
      geom_point(aes(x = X1, y = X2,size = !!sym(size),fill = !!sym(fill)),pch = 21, data =  plotnode) + scale_colour_brewer(palette = "Set1") +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
      labs( title = paste(layout,"network",sep = "_")) + theme_void()
    pnet

    if (label == TRUE ) {
      pnet <- pnet +  geom_text(aes(X1, X2,label= !!sym(lab)), data = plotnode)

    }

    plotname = paste(path,"/network",layout,".pdf",sep = "")
    ggsave(plotname, pnet, width = width, height =height)

    plots[[aa]] = pnet
    #----
    rand.g<- erdos.renyi.game(length(V(igraph)), length(E(igraph)),type = c("gnm"))
    # degree_distribution

    ## compared the degree of random network and this network
    data1 = data.frame(network= degree_distribution(igraph, cumulative = FALSE),group = "E–R network",ID = c(1:length(degree_distribution(igraph, cumulative = FALSE))))
    data2 = data.frame(network = degree_distribution(rand.g, cumulative = FALSE) ,group = "network",ID = c(1:length(degree_distribution(rand.g, cumulative = FALSE) )))
    data = rbind(data1,data2)
    p1 <- ggplot(data) +geom_point(aes(x = ID,y = network,group =group,fill = group),pch = 21,size = 2) +
      geom_smooth(aes(x = ID,y = network,group =group,color = group))+
      theme_bw()
    plotname = paste(path,"/Power_law_distribution_",layout,".pdf",sep = "")
    ggsave(plotname, p1, width = width, height =height)

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
  p  = ggarrange(plotlist = plots, common.legend = TRUE, legend="right",ncol = ncol,nrow = nrow)
  return(list(p,y))

}



