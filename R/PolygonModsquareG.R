#' Construct a network layout. Arrange network nodes to different locations according to grouping
#'
#' @param cor Correlation matrix
#' @param nodeGroup Classification information of network nodes.Group according to actual requirements, see example
#' @param r1 big cluster r
#' @param N break of each cluster
#' @param cut line need to cut off
#' @examples
#' data
#' data(ps)
#' result = corMicro (ps = ps,N = 0.02,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # building the node group
#' netClu = data.frame(ID = row.names(cor),group =rep(1:3,length(row.names(cor)))[1:length(row.names(cor))] )
#' netClu$group = as.factor(netClu$group)
#' result2 = PolygonModsquareG(cor = cor,nodeGroup =netClu,r1 = 1,N = 1.1,cut = 3,line = FALSE)
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




PolygonModsquareG <- function(cor = cor,nodeGroup =netClu,r1 = 1,N = 1.1,cut = 3){

  mod = as.data.frame(table(nodeGroup$group)) %>%
    arrange(desc(Freq))
  r = c(r1,rep(0,(dim(mod)[1]-1)))
  for (i in 2:dim(mod)[1]) {
    r[i] = r1*mod[i,2]/mod[1,2]

  }
  #
  # r * 2 + 1

  data = data.frame(x = c((cumsum(r * 2 + N))[1]  ,(cumsum(r * 2 + N))[-1]),y = 0 )
  data
  if (cut != 1) {
    # A = r * 2 + N
    num = 1
    # i = 2
    A = r * N
    # B = data$x
    r * 2 + N
    for (i in num:length(A)) {
      if (sum(A[num:i]) > sum((r * N))/cut | i == length(A) ) {
        data$y[num:i] = min(data$x[num:i])
        data$x[num:i] =  data$x[num:i] -  data$x[num]
        print(i)
        num = i + 1




      }
    }

  }
  data
  da = data

  nodeGroup$group = factor(nodeGroup$group ,levels = mod$Var1)



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

    as$ID
    as$ID = as.character( as$ID)

    # Calculation of a single circular coordinate

    if (length(as$ID)!=1 ) {
      m = cor[as$ID,as$ID]

      d  =m
      d <- as.edgelist.sna(d)
      # if (is.list(d))
      # d <- d[[1]]
      n <- attr(d, "n")
      #
      s = r[i]


      data = cbind(sin(2 * pi * ((0:(n - 1))/n))*s +da[i,1], cos(2 * pi * ((0:(n - 1))/n))*s +da[i,2])

      data =as.data.frame(data)
      row.names(data ) = row.names(m)
      data$elements = row.names(data )
      colnames(data)[1:2] = c("X1","X2")

    }



    head(data)

    # ggplot(data) + geom_point(aes(x = X1,y = X2))

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

