#' Construct a network layout. Randomly distributed module positions
#'
#' @param cor Correlation matrix
#' @param nodeGroup Classification information of network nodes
#' @examples
#' data
#' data(ps)
#' result = corMicro (ps = ps,N = 0.02,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # Extract tax table for grouping
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
#' # Calculate the layout of netwok, extract the coordinate of node
#' result2 = randomClusterG (cor = cor,nodeGroup =netClu )
#' node = result2[[1]]
#'
#'
#' @return result2 Which contains 2 data.frame. Result2[[1]], consists of OTU and its corresponding coordinates.
#' result2[[2]], consists of the network center coordinates of each group
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export



randomClusterG = function(cor = cor,nodeGroup =netClu ){

  num = length(levels(nodeGroup$group))
  xs = as.data.frame(table(nodeGroup$group))
  r = xs$Freq/10

  #--Calculate the total edge radius and#-----
  rtotal = sum(r)

  #---Start the coordinates of the clustering module#-----
  whichnum = dim(combn(num,2))[2]
  # Pairwise comparison, the radius is greater than a certain threshold is TRUE, all TRUE is qualified#---------
  A = rep(FALSE,whichnum)
  A

  while (all(A) == FALSE) {



    #--Randomly generate a set of coordinates-these coordinate ranges are limited to the sum of all cluster radii#----
    da = data.frame(x =sample((0:rtotal),num ),y =  sample((rtotal):0,num ))
    head(da)


    for (i in 1:whichnum) {
      #--Comparison between the two groups
      cs = combn(num,2)
      cs
      # i =13
      cs[,i]
      # Extract the coordinates of the two sets of cluster center coordinates
      x1 = da[cs[,i],][1,1]
      x2 = da[cs[,i],][2,1]
      y1 = da[cs[,i],][1,2]
      y2 = da[cs[,i],][2,2]
      # Determine whether the distance between these two coordinates is greater than the sum of the two clustering radii
      AS= (x1-x2)^2 +(y1-y2)^2
      A[i]= sqrt(AS)> sum(r[cs[,i]])

    }

  }

  #-Start calculating layout#-----
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


    #-----------Calculation of a single circular coordinate

    if (length(as$ID)!=1 ) {
      m = cor[as$ID,as$ID]

      d  =m
      d <- as.edgelist.sna(d)
      # if (is.list(d))
      # d <- d[[1]]
      n <- attr(d, "n")
      # 提取半径
      s = r[i]



      # if (i == 1) {
      # sori = s*2
      data = cbind(sin(2 * pi * ((0:(n - 1))/n))*s +da[i,1], cos(2 * pi * ((0:(n - 1))/n))*s +da[i,2])

      # }
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
