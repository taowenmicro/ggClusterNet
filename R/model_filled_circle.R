

model_filled_circle <- function(cor = cor,
                                culxy =TRUE,
                                da = NULL,# 数据框，包含x,和y列
                                nodeGroup =netClu,
                                seed = 10,
                                mi.size = 0.5,
                                zoom = 1
                                ){

  if (culxy ==TRUE) {
    #--多边形布局的核心位置的计算#-------------
    num = length(levels(nodeGroup$group))
    #---Extract the number of nodes in each group and define the circle radius according to the number
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
    da = data.frame(x = x,y = y)
  } else{
    da = da
  }

  # ggplot(da) + geom_point(aes(x,y))
  #--计算每个部分的排布方式#----------
  nodeGroup$group = as.factor(nodeGroup$group)
  #-Start calculating layout
  # i= 5
  for (i in 1:length(levels(nodeGroup$group))) {

    #--Extract all otu in this group
    as = dplyr::filter(nodeGroup, group == levels(nodeGroup$group)[i])
    as$ID = as.character( as$ID)
    #-如果这个模块只有一个，那就放到圆心位置
    if (length(as$ID) == 1) {
      data = cbind(da[i,1],da[i,2] )
      data =as.data.frame(data)
      row.names(data ) = as$ID
      data$elements = row.names(data)
      colnames(data)[1:2] = c("X1","X2")
    }



    #Calculation of a single circular coordinate
    if (length(as$ID)!= 1 ) {

      m = cor[as$ID,as$ID]
      num.node <- dim(m)[1]

      # num.node  = 6
      if (num.node <=20&num.node >=7) {
        n = 1
      } else if(num.node >1 &num.node <7){
        n = 0
      } else {
        for (N in 1: num.node) {
          A = 1 + (7*(N + 1)*N )/2 - N
          if (A >= num.node) {
            break
          }
          n = N - 1
          print(n)
        }

      }

      if (n >= 1) {
        # n = (sqrt((num.node-1)/3) - 1) %>% floor()
        wai.mode = num.node - (1 + (7*(n + 1)*n )/2 - n)
        dat = data.frame(x = 0,y = 0)

        for (ii in 1:n) {
          t <- seq(0, 2*pi, length.out = 7*ii)
          t = t[-1]
          x <- sin(t)*ii
          y <- cos(t)*ii
          add = data.frame(x = x,y = y)
          dat = rbind(dat,add)

          if (ii== n) {
            ii = ii + 1
            t <- seq(0, 2*pi, length.out = (wai.mode + 1))
            t = t[-1]
            x <- sin(t)*ii
            y <- cos(t)*ii
            add = data.frame(x = x,y = y)
            dat = rbind(dat,add)
          }

        }
      } else if(n == 0){
        dat = data.frame(x = 0,y = 0)
        t <- seq(0, 2*pi, length.out = (num.node))
        t = t[-1]
        x <- sin(t)*mi.size
        y <- cos(t)*mi.size
        add = data.frame(x = x,y = y)
        dat = rbind(dat,add)



      }


      dim(dat)
      # ggplot(dat) + geom_point(aes(x,y))
      data =data.frame(x = dat$x + da[i,1],y = dat$y + da[i,2])
      # ggplot(data) + geom_point(aes(x,y))
      row.names(data) = row.names(m)
      data$elements = row.names(m)
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
  da$group = levels(nodeGroup$group)
  row.names(da) = levels(nodeGroup$group)
  return(list(plotcord,da))

}
# ggplot(oridata) + geom_point(aes(X1,X2))


