#' Use the layout in the sne package to calculate the visual layout of the network
#'
#' @param cor Correlation matrix
#' @param layout select layout to building sub cluster，could br on of the sna package layout, such as "circle","adj","circrand","eigen","random".
#' @param nodeGroup group you must imput
#' @param zoom Set the distance between modules
#' @examples
#' data
#' data(ps)
#' result = corMicro (ps = ps,N = 0.02,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # building the node group
#' netClu = data.frame(ID = row.names(cor),group =rep(1:3,length(row.names(cor)))[1:length(row.names(cor))] )
#' netClu$group = as.factor(netClu$group)
#' result2 = ranSNEClusterG (cor=  cor,layout ="circle",nodeGroup = netClu)
#' node = result2[[1]]
#'
#'
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


ranSNEClusterG = function(cor=  cor,layout ="circle",nodeGroup = netClu,zoom = 1){

  num = length(levels(nodeGroup$group))
  xs = as.data.frame(table(nodeGroup$group))
  r = rep(1,length(xs$Freq))
  rtotal = sum(r) *zoom

  whichnum = dim(combn(num,2))[2]
  # whichnum = num
  A = rep(FALSE,whichnum)
  A
  while (all(A) == FALSE) {
    da = data.frame(x =sample((0:rtotal),num ),y =  sample((rtotal):0,num ))
    for (i in 1:whichnum) {

      cs = combn(num,2)
      cs
      # i =13
      cs[,i]
      #
      x1 = da[cs[,i],][1,1]
      x2 = da[cs[,i],][2,1]
      y1 = da[cs[,i],][1,2]
      y2 = da[cs[,i],][2,2]
      #
      AS= (x1-x2)^2 +(y1-y2)^2
      A[i]= sqrt(AS)> sum(r[cs[,i]])
    }
  }


  #
  for (i in 1:length(levels(nodeGroup$group))) {

    #-
    as = dplyr::filter(nodeGroup, group == levels(nodeGroup$group)[i])


    if (length(as$ID) == 1) {
      data = cbind(da[i,1],da[i,2] )

      data =as.data.frame(data)
      row.names(data ) = as$ID
      data$elements = row.names(data )
      colnames(data)[1:2] = c("X1","X2")
    }

    as$ID = as.character( as$ID)


    #-

    if (length(as$ID)!=1 ) {

      m = cor[as$ID,as$ID]
      #-------尝试全部的layout#----
      if (layout == "circle") {
        data<- data.frame(sna::gplot.layout.circle(m, NULL))
      }

      if (layout == "kamadakawai") {
        data <- data.frame(gplot.layout.adj(m, NULL))
      }
      if (layout == "adj") {
        data <- data.frame(gplot.layout.kamadakawai(m, NULL))
      }
      if (layout == "circrand") {
        data <- data.frame(gplot.layout.circrand(m, NULL))
      }
      if (layout == "eigen") {
        data <- data.frame(gplot.layout.eigen(m, NULL))
      }
      if (layout == "geodist") {
        data <- data.frame(gplot.layout.geodist(m, NULL))
      }
      if (layout == "hall") {
        data <- data.frame(gplot.layout.hall(m, NULL))
      }

      if (layout == "mds") {
        data <- data.frame(gplot.layout.mds(m, NULL))
      }

      if (layout == "princoord") {
        data <- data.frame(gplot.layout.princoord(m, NULL))
      }
      if (layout == "random") {
        data <- data.frame(gplot.layout.random(m, NULL))
      }
      if (layout == "rmds") {
        data <- data.frame(gplot.layout.rmds(m, NULL))
      }
      if (layout == "segeo") {
        data <- data.frame(gplot.layout.segeo(m, NULL))
      }
      if (layout == "seham") {
        data <- data.frame(gplot.layout.seham(m, NULL))
      }
      if (layout == "spring") {
        data <- data.frame(gplot.layout.spring(m, NULL))
      }
      if (layout == "springrepulse") {
        data <- data.frame(gplot.layout.springrepulse(m, NULL))
      }
      if (layout == "target") {
        data <- data.frame(gplot.layout.target(m, NULL))
      }


      head(data)

      #------After the calculation is completed, remember that this is calculated according to the origin coordinate of 0, we need to modify it to our coordinates
      data$X1 = data$X1 +da[i,1]
      data$X2 = data$X2 +da[i,2]
      data$elements <- colnames(m)
      row.names(data ) = data$elements

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



