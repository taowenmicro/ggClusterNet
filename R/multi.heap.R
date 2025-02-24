multi.heap = function(
    zoom.x = 20,
    zoom.y = 10,
    matrix.env = list(env1 = env1,
                      env2 = env1,
                      env3 = env1,
                      env4 = env1
    )




){
  for (i in 1:length(names(matrix.env))) {
    tem = matrix.env[[names(matrix.env)[i]]]
    res1  = cor_heat(data = tem,
                     text = "down",
                     label_text = F,
                     # label = c(3,4),
                     show_labels = c("left",  "right"),# "right", "top",
                     label_rotation = list(left = 0, right = 0, top = 0, bottom = 0), # 文字角度
                     label_padding = 0.5 )
    tem2 = res1[[2]]
    tem2$Class = names(matrix.env)[i]

    if (i ==1) {
      tem0 = tem2
    } else{
      tem0 = rbind(tem0,tem2)
    }


  }


  if (length(names(matrix.env)) == 2) {
    j = 2
    num = tem0$Class %>% unique()
    tem0$x[tem0$Class %in% num[j]]  = tem0$x[tem0$Class %in% num[j]]  + zoom.x
    tem0$y[tem0$Class %in% num[j]]  = tem0$y[tem0$Class %in% num[j]]  + zoom.y

  }else if (length(names(matrix.env)) == 3){
    j = 2
    num = tem0$Class %>% unique()
    tem0$x[tem0$Class %in% num[j]]  = tem0$x[tem0$Class %in% num[j]]  + zoom.x
    tem0$y[tem0$Class %in% num[j]]  = tem0$y[tem0$Class %in% num[j]]  + 0
    j = 3
    num = tem0$Class %>% unique()
    tem0$x[tem0$Class %in% num[j]]  = tem0$x[tem0$Class %in% num[j]]  + zoom.x/2
    tem0$y[tem0$Class %in% num[j]]  = tem0$y[tem0$Class %in% num[j]]  + zoom.y



  }else if (length(names(matrix.env)) == 4){
    j = 2
    num = tem0$Class %>% unique()
    tem0$x[tem0$Class %in% num[j]]  = tem0$x[tem0$Class %in% num[j]]  + zoom.x
    tem0$y[tem0$Class %in% num[j]]  = tem0$y[tem0$Class %in% num[j]]  + 0
    j = 3
    num = tem0$Class %>% unique()
    tem0$x[tem0$Class %in% num[j]]  = tem0$x[tem0$Class %in% num[j]]  + zoom.x/2
    tem0$y[tem0$Class %in% num[j]]  = tem0$y[tem0$Class %in% num[j]]  + zoom.y
    j = 4
    num = tem0$Class %>% unique()
    tem0$x[tem0$Class %in% num[j]]  = tem0$x[tem0$Class %in% num[j]]  + zoom.x/2
    tem0$y[tem0$Class %in% num[j]]  = tem0$y[tem0$Class %in% num[j]]  + -zoom.y

  }



  gradient_colors <- c("#2E8B57", "white", "#CD5C5C")
  gradient_breaks <- seq(-1, 1, 0.5)
  p =  ggplot( data = tem0, aes(x = x, y = y)) +
    # 全矩阵背景网格
    geom_tile(fill = "white", color = "gray90", linewidth = 0.3)+
    geom_tile(

      aes(fill = cor, width = abs(cor)*0.9, height = abs(cor)*0.9),
      color = "gray90",
      linewidth = 0.3) + theme_void() +
    # 颜色映射设置
    scale_fill_gradientn(
      colors = gradient_colors,
      values = scales::rescale(gradient_breaks),
      limits = c(-1, 1),
      breaks = gradient_breaks) +

    theme(plot.margin=unit(rep(2,4),'cm'))  +
    # 图例设置
    guides(fill = guide_colorbar(
      title = "Correlation",
      barwidth = 1,
      barheight = 10,
      frame.colour = "black",
      ticks.colour = "black"))

  return(list(p,tem0))
}


