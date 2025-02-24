# # remotes::install_github("taowenmicro/ggClusterNet")
# library(ggClusterNet)
# library(reshape2)
#
# otu1 =   as.data.frame(t(ggClusterNet::vegan_otu(ps)))
# otu2 =   as.data.frame(t(ggClusterNet::vegan_otu(ps)))
# otu3 =   as.data.frame(t(ggClusterNet::vegan_otu(ps)))
#
# tabOTU = list(bac = otu1,
#               bac2 =otu2,
#               bac3= otu3)
#
# matrix.line = list(right = c(names(tabOTU)[1:2]),
#                    left = c(names(tabOTU)[1:2]),
#                    up = c(names(tabOTU)[2:3]),
#                    bottom = c(names(tabOTU)[2]))
#
# rep = MetalTast (env.dat = env1, tabOTU = tabOTU)
# repR = rep[c(-seq(from=1,to=dim(rep)[2],by=2)[-1])]
# repP = rep[seq(from=1,to=dim(rep)[2],by=2)]
#
# mantel = cbind(repR, repP)
# mantel = mantel[, !duplicated(colnames(mantel))]
#
# res  = cor_heat(data = env1,
#                 text = "down",
#                 label_text = F,
#                 # label = c(3,4),
#                 show_labels = c("left",  "right"),# "right", "top",
#                 label_rotation = list(left = 0, right = 0, top = 0, bottom = 0), # 文字角度
#                 label_padding = 0.5 )
# p = res[[1]]
# p
#
# dat = res[[2]]
# head(dat)

cor_heat = function(data = env1,
                    text = "down",
                    label_text =  FALSE,
                    # label = c(3,4),
                    show_labels = c("left",  "bottom"),# "right", "top",

                    label_rotation = list(left = 0, right = 0, top = 0, bottom = 0), # 文字角度
                    label_padding = 0.5

){

  occor = psych::corr.test(data , use = "pairwise",adjust = "fdr")
  occor.r = occor$r
  occor.p = occor$p

  cor_matrix <- as.matrix(occor.r)

  p_matrix <- as.matrix(occor.p)


  if(label_text ==  FALSE){
    # 步骤2：创建整数坐标系统
    n_vars <- ncol(cor_matrix)

    coord_map <- data.frame(
      var = colnames(cor_matrix),
      x = 1:n_vars,
      y = n_vars:1  # y轴逆序排列
    )

    # 步骤3：合并数据
    plot_data <- reshape2::melt(cor_matrix, value.name = "cor") %>%
      left_join(coord_map %>% select(var, y_coord = y),
                by = c("Var1" = "var")) %>%
      left_join(coord_map %>% select(var, x_coord = x),
                by = c("Var2" = "var")) %>%
      select(Var1, Var2, x = x_coord, y = y_coord, cor)

    gradient_colors <- c("#2E8B57", "white", "#CD5C5C")
    gradient_breaks <- seq(-1, 1, 0.5)

    p = ggplot(plot_data, aes(x = x, y = y)) +
      # 全矩阵背景网格
      geom_tile(fill = "white", color = "gray90", linewidth = 0.3)+
      geom_tile(
        data = plot_data,
        aes(fill = cor, width = abs(cor)*0.9, height = abs(cor)*0.9),
        color = "gray90",
        linewidth = 0.3) + theme_void()

    p = p+
      # 颜色映射设置
      scale_fill_gradientn(
        colors = gradient_colors,
        values = scales::rescale(gradient_breaks),
        limits = c(-1, 1),
        breaks = gradient_breaks) +
      theme(
        axis.title = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA)) +
      theme(plot.margin=unit(rep(2,4),'cm'))+
      # 图例设置
      guides(fill = guide_colorbar(
        title = "Correlation",
        barwidth = 1,
        barheight = 10,
        frame.colour = "black",
        ticks.colour = "black"))
  }
  else if (label_text ==  TRUE){

    # 生成坐标映射
    n_vars <- ncol(cor_matrix)
    coord_map <- data.frame(
      var = colnames(cor_matrix),
      x = 1:n_vars,           # x轴正序排列
      y = 1:n_vars            # y轴正序排列（后续通过scale_y_reverse反转显示）
    )


    plot_data <- reshape2::melt(cor_matrix, value.name = "cor") %>%
      left_join(melt(p_matrix, value.name = "p"), by = c("Var1", "Var2"))  %>%
      left_join(coord_map %>% select(var, y_coord = y),
                by = c("Var1" = "var")) %>%
      left_join(coord_map %>% select(var, x_coord = x),
                by = c("Var2" = "var")) %>%
      select(Var1, Var2, x = x_coord, y = y_coord, cor,p)  %>%
      # 排序保证坐标系正确
      #  arrange(x, desc(y)) %>%
      mutate(
        cell_type = case_when(
          x < y ~ "upper",
          x > y ~ "lower",
          TRUE ~ "diag"
        ),
        sig = case_when(
          p < 0.001 ~ "***",
          p < 0.01 ~ "**",
          p < 0.05 ~ "*",
          TRUE ~ ""
        ),
        label = sprintf("%.2f%s", cor, sig)
      )

    # 自定义颜色梯度
    gradient_colors <- c("#2E8B57", "white", "#CD5C5C")
    gradient_breaks <- seq(-1, 1, 0.5)
    env_levels <- as.factor(colnames(env1))


    # 生成热图
    if (text == "down"){

      p = ggplot(plot_data, aes(x = x, y =y)) +
        # 全矩阵背景网格
        geom_tile(fill = "white", color = "gray90", linewidth = 0.3) +
        # 上三角热图区块
        geom_tile(
          data = filter(plot_data, cell_type == "upper"),
          aes(fill = cor, width = abs(cor)*0.9, height = abs(cor)*0.9),
          color = "gray90",
          linewidth = 0.3
        ) +
        theme_void()+
        # 下三角数值标签
        geom_text(
          data = filter(plot_data, cell_type == "lower"),
          aes(label = label),
          color = "black",
          size = 3.5
        ) +
        # 对角线处理
        geom_tile(
          data = filter(plot_data, cell_type == "diag"),
          fill = "gray95",
          color = "gray90"
        ) +
        # 颜色映射设置
        scale_fill_gradientn(
          colors = gradient_colors,
          values = scales::rescale(gradient_breaks),
          limits = c(-1, 1),
          breaks = gradient_breaks
        ) +
        # 坐标轴优化
        coord_equal() +
        # scale_x_discrete(position = "top") +
        # scale_y_discrete(limits = rev) +
        # 主题美化
        theme_minimal(base_size = 12) +
        theme(
          axis.title = element_blank(),
          #  axis.text.x = element_text(angle = 45, hjust = 0, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA)
        ) +
        theme(plot.margin=unit(rep(2,4),'cm'))+
        # 图例设置
        guides(fill = guide_colorbar(
          title = "Correlation",
          barwidth = 1,
          barheight = 10,
          frame.colour = "black",
          ticks.colour = "black"
        ))
      p

    }
    else if (text == "up") {

      p = ggplot(plot_data, aes(x = x, y =y)) +
        # 全矩阵背景网格
        geom_tile(fill = "white", color = "gray90", linewidth = 0.3) +
        # 上三角文字区块
        geom_tile(
          data = filter(plot_data, cell_type == "lower"),
          aes(fill = cor, width = abs(cor)*0.9, height = abs(cor)*0.9),
          color = "gray90",
          linewidth = 0.3
        ) +
        theme_void()+
        # 下三角数值标签
        geom_text(
          data = filter(plot_data, cell_type == "upper"),
          aes(label = label),
          color = "black",
          size = 3.5
        ) +
        # 对角线处理
        geom_tile(
          data = filter(plot_data, cell_type == "diag"),
          fill = "gray95",
          color = "gray90"
        ) +
        # 颜色映射设置
        scale_fill_gradientn(
          colors = gradient_colors,
          values = scales::rescale(gradient_breaks),
          limits = c(-1, 1),
          breaks = gradient_breaks
        ) +
        #   # 坐标轴优化
        coord_equal() +
        #   scale_x_discrete(position = "top") +
        #   scale_y_discrete(limits = rev) +
        # 主题美化
        theme_minimal(base_size = 12) +
        theme(
          axis.title = element_blank(),
          # axis.text.x = element_text(angle = 45, hjust = 0, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA)
        ) +
        # 图例设置
        guides(fill = guide_colorbar(
          title = "Correlation",
          barwidth = 1,
          barheight = 10,
          frame.colour = "black",
          ticks.colour = "black"
        ))
      p
    }
  }



  p =  p + theme(axis.text.y = element_blank(),axis.text.x = element_blank())

  # 生成标签定位数据
  label_pos <- list(
    left = data.frame(
      x = 1 - label_padding,
      y = coord_map$y,
      label = coord_map$var,
      hjust = 1,
      vjust = 0.5
    ),
    right = data.frame(
      x = n_vars + label_padding,
      y = coord_map$y,
      label = coord_map$var,
      hjust = 0.5,
      vjust = 0.5
    ),
    top = data.frame(
      x = coord_map$x,
      y = n_vars + label_padding,
      label = coord_map$var,
      hjust = 0.5,
      vjust = 0
    ),
    bottom = data.frame(
      x = coord_map$x,
      y = 1 - label_padding,
      label = coord_map$var,
      hjust = 0.5,
      vjust = 1.5
    )
  )

  # 选择需要显示的标签
  labels_to_show <- do.call(rbind, lapply(show_labels, function(pos) {
    df <- label_pos[[pos]]
    df$position <- pos
    df$angle <- label_rotation[[pos]]
    df
  }))

  # 添加标签
  if(nrow(labels_to_show) > 0) {
    p1 <- p + geom_text(
      data = labels_to_show,
      aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
      angle = labels_to_show$angle,
      size = 3.5 )
  }

  p1 =p1 +
    coord_cartesian(
      clip = "off",          # 允许标签溢出绘图区
      expand = FALSE         # 完全禁用扩展
    )

  return(list(
    plot = p1,
    plot_data = plot_data
  ))

}

