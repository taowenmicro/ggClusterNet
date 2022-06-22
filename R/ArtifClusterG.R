#' Manually specify coordinates to display the network
#'
#' @param cor Correlation matrix
#' @param nodeGroup Classification information of network nodes.Group according to actual requirements, see example
#' @param r Radius of each submodule
#' @param da The coordinates of each submodule
#' @examples
#' data(ps)
#' result = corMicro (ps = ps,N = 100,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # Extract taxonomy table
#' ps_net = result[[3]]
#' # set group of experimental
#' netClu = data.frame(ID = row.names(tax_table),group =rep(1:5,length(row.names(tax_table)))[1:length(row.names(tax_table))] )
#' netClu$group = as.factor(netClu$group)
#' tax_table = ps_net %>% vegan_tax() %>%
#' as.data.frame()
#' xs = as.data.frame(table(netClu$group))
#' # set the radicus
#' r = rep(4,length(xs$Freq))
#' # Set the coordinates manually
#' ax1 = c(120,0)
#' ax2 = c(130,-30)
#' ax3 = c(140,-70)
#' ax4 = c(130,-110)
#' ax5 = c(120,-140)
#' da = rbind(ax1,ax2,ax3,ax4,ax5)
#' # Calculate network layout
#' result2 = ArtifCluster(cor = cor,nodeGroup =netClu,r = r,da =da)
#' node = result2[[1]]
#' @return list which contains node position coordinates
#' @author Contact: Tao Wen \email{taowen@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn} yongxin liu \email{yxliu@@genetics.ac.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Tao Wen#, Penghao Xie#, Shengdie Yang, Guoqing Niu, Xiaoyu Liu, Zhexu Ding, Chao Xue, Yong-Xin Liu *, Qirong Shen, Jun Yuan*
#' ggClusterNet: an R package for microbiome network analysis and modularity-based multiple network layouts
#' iMeta 2022,DOI: \url{doi: 10.1002/imt2.32}
#' @export



ArtifCluster = function(cor = cor,nodeGroup =netClu,r = r,da =da){
  nodeGroup$group = as.factor(nodeGroup$group)
  for (i in 1:length(levels(nodeGroup$group))) {
    #--Extract all otu in this group
    as = dplyr::filter(nodeGroup, group == levels(nodeGroup$group)[i])
    if (length(as$ID) == 1) {
      data = cbind(da[i,1],da[i,2] )
      data =as.data.frame(data)
      row.names(data ) = as$ID
      data$elements = row.names(data )
      colnames(data)[1:2] = c("X1","X2")
    }
    as$ID = as.character( as$ID)
    # Calculation of a single circular coordinate
    if (length(as$ID)!=1 ) {
      m = cor[as$ID,as$ID]
      d  =m
      d <- as.edgelist.sna(d)
      n <- attr(d, "n")
      s = r[i]
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
