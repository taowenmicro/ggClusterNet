#' Construct a network layout. Arrange network nodes to different locations according to grouping
#'
#' @param cor Correlation matrix
#' @param nodeGroup Classification information of network nodes.Group according to actual requirements, see example
#' @param order node order, TURE or FALSE
#' @examples
#' data
#' data(ps)
#' result = corMicro (ps = ps,N = 0.02,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # building the node group
#' netClu = data.frame(ID = row.names(cor),group =rep(1:3,length(row.names(cor)))[1:length(row.names(cor))] )
#' netClu$group = as.factor(netClu$group)
#' result2 = PolyRdmNotdCirG(cor = cor,nodeGroup =  netClu )
#'
#'
#' @return result2 Which contains 2 lists.Result2[[1]], consists of OTU and its corresponding coordinates.
#' Result2[[2]], consists of the network center coordinates of each group
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


PolyRdmNotdCirG <- function(cor = cor,nodeGroup = netClu, order = FALSE){

  num = length(levels(nodeGroup$group))
  #---Extract the number of nodes in each group and define the circle radius according to the number
  xs = as.data.frame(table(nodeGroup$group))
  r = xs$Freq/10
  #--Calculate angle according to group
  arg = seq(0,360,360/(length(r)))
  i = 1
  rsum = sum(r)/2
  x= rep(0,length(r))
  y = rep(0,length(r))

  for (i in 1:length(r)) {

    x[i] = (rsum)* sin(arg[i]* 3.14/180)

    y[i] = (rsum)* cos(arg[i]* 3.14/180)
  }

  da = data.frame(x = x,y = y)

  nodeGroup$group = as.factor(nodeGroup$group)
  #-Start calculating layout#-----
  for (i in 1:length(levels(nodeGroup$group))) {

    #--Extract all otu in this group
    as = dplyr::filter(nodeGroup, group == levels(nodeGroup$group)[i])
    if (length(as$ID) == 1) {
      data = cbind(da[i,1],da[i,2] )
      data =as.data.frame(data)
      row.names(data ) = as$ID
      data$elements = row.names(data)
      colnames(data)[1:2] = c("X1","X2")
    }
    as$ID = as.character( as$ID)

    #Calculation of a single circular coordinate
    if (length(as$ID)!=1 ) {
      m = cor[as$ID,as$ID]

      #
      # d  =m
      # d <- as.edgelist.sna(d)
      # n <- attr(d, "n")
      # # length(colnames(m))
      # # Extract radius
      # # s = r[i]
      # s = mean(r)
      # data = cbind(sin(2 * pi * ((0:(n - 1))/n))*s +da[i,1], cos(2 * pi * ((0:(n - 1))/n))*s +da[i,2])

      if (order== TRUE) {
        packing <- packcircles::circleProgressiveLayout(rep(1,length(colnames(m))))

      } else {
        packing <- packcircles::circleProgressiveLayout(runif(n=length(colnames(m))))

      }




      data <- packcircles::circleLayoutVertices(packing)  %>% dplyr::group_by(id) %>%
        dplyr::summarise(x = mean(x),y = mean(y))  %>%
        dplyr::select(-id)  %>%
        as.data.frame()
     # head(data)
       data <- data.frame(x = data$x + da[i,1],y = data$y + da[i,2])


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
