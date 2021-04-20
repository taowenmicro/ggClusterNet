#' Calculate a bipartite network for microbiome data
#'
#' @param table ven output
#' @param distance core distance to group
#' @param distance2
#' @param distance3
#' @param order TRUE or FEASE
#'
#' @param flour Is the arrangement of points discrete? TRUE or FEASE would be selected
#' @examples
#' data(ps)
#' ps_sub = filter_taxa(ps, function(x) sum(x ) > 20 , TRUE)
#' ps_sub = filter_taxa(ps_sub, function(x) sum(x ) < 30 , TRUE)
#' ps_sub
#'
#' result = div_network(ps_sub,num = 6)
#'
#' edge = result[[1]]
#' head(edge)
#' result <- div_culculate(table = result[[3]],distance = 1.1,distance2 = 1.5,distance3 = 1.3,order = FALSE)
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


div_culculate <- function(table = result[[3]],distance = 1.1,distance2 = 1.5,distance3 = 1.3,order = FALSE){

  all = table[rowSums(table) == dim(table)[2],]
  allnum = dim(table[rowSums(table) == dim(table)[2],])[1]
  N = allnum
  if (N == 0) {
    print("all cover was 0,so,can not continue")
  }else {
    if (order== TRUE) {
      packing <- packcircles::circleProgressiveLayout(rep(1,N))

    }
    if (order== FALSE) {
      packing <- packcircles::circleProgressiveLayout(runif( min = 1, max = 10,n=N))

    }

    data <- packcircles::circleLayoutVertices(packing)  %>% dplyr::group_by(id) %>%
      dplyr::summarise(x = mean(x),y = mean(y))  %>%
      dplyr::select(-id)  %>%
      as.data.frame()
    dim(data)
    r0 = max(data$x) - min(data$x)
    row.names(data ) = row.names(all)
    data$elements = row.names(data )
    colnames(data)[1:2] = c("X1","X2")
    allxy = data


    #---culculate da
    r = distance*r0
    #--Calculate angle according to group
    arg = seq(0,360,360/(dim(table)[2]))
    x= rep(0,dim(table)[2])
    y = rep(0,dim(table)[2])
    for (i in 1:dim(table)[2]) {

      x[i] = r* sin(arg[i]* 3.14/180)

      y[i] = r* cos(arg[i]* 3.14/180)
    }

    da0 = data.frame(x = x,y = y)

    #---------
    r = distance2*r0
    x= rep(0,dim(table)[2])
    y = rep(0,dim(table)[2])
    for (i in 1:dim(table)[2]) {
      x[i] = r* sin(arg[i]* 3.14/180)
      y[i] = r* cos(arg[i]* 3.14/180)
    }
    da1 = data.frame(x = x,y = y)
    da1
    # i = 1
    for (i in 1:dim(table)[2]) {

      table = table[rowSums(table) != dim(table)[2],]
      N = length(table[,i][table[,i] != 0])

      if (N != 0) {
        if (order== TRUE) {
          packing <- packcircles::circleProgressiveLayout(rep(1,N))
        } else {
          packing <- packcircles::circleProgressiveLayout(runif(min = 1, max = 10,n=N))
        }
        data <- packcircles::circleLayoutVertices(packing)  %>% dplyr::group_by(id) %>%
          dplyr::summarise(x = mean(x),y = mean(y))  %>%
          dplyr::select(-id)  %>%
          as.data.frame()
        data <- data.frame(x = data$x + da1[i,1] * distance3,y = data$y + da1[i,2]* distance3)
        row.names(data ) = row.names(table[table[,i] != 0,])
        data$elements = row.names(data )
        colnames(data)[1:2] = c("X1","X2")
      } else{
        data = NULL
      }

      if (i == 1) {
        oridata = data
      }
      if (i != 1) {
        oridata = rbind(oridata,data)
      }
    }

    oridata = rbind(oridata,allxy)
    head(oridata )

    rownames(da0) = colnames(table)
    da0$elements = row.names(da0 )
    colnames(da0)[1:2] = c("X1","X2")
    # plotdata = rbind(oridata,da0)
    plotdata=  oridata

    edg <-plotdata %>%
      dplyr::rename(X1 = X1,Y1 = X2,OTU = elements) %>%
      # select(X) %>%
      dplyr::right_join(edge,by = c("OTU"="source"))

    head(edg)
    head(plotdata)
    edge2 <-da0 %>%
      dplyr::rename(X2 = X1,Y2 = X2,Group = elements) %>%
      dplyr::right_join(edg,by = c("Group"="target"))
    return(list(edge2,plotdata,da0))
  }


}

