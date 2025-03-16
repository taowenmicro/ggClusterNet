
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
#
# res2 = collinkh (zoom  = 5,# 微生物数据和热图距离
#                  offset1=1,  # 微调连线和热图距离
#                  corva = -0.2 ,
#                  angle = 40,
#                  sig = F,
#                  p.thur = 0.3,
#                  heat_coords = dat,
#                  mantel = mantel,
#                  tabOTU = tabOTU,
#                  p = res[[1]],
#                  matrix.line =  matrix.line
# )
#
#
#
# res2[[1]]


collinkh2 = function(zoom  = 4,# 微生物数据和热图距离
                    offset1=0.7,  # 微调连线和热图距离
                    corva = 0.2 ,
                    angle = 40,
                    sig = F,
                    heat_coords = heat_coords,
                    mantel = mantel,
                    matrix.line =  matrix.line,
                    p = plotggplot,
                    p.thur = 0.05
){

  if (sig) {
    mantel_p  = "sig"
  } else{
    mantel_p = "nosig"
  }
  input_p = p

  heat_coords=heat_coords %>%
    dplyr::rename(
      x_env = x,
      y_env = y
      #  env_var = Var1  # 可选: 更清晰的列名
    )

  all_dfs <- list()

  id.right = matrix.line$right
  if(!is.null(id.right)) {

    microbe_y_range <- range(heat_coords$y_env)
    # tabOTU_l = tabOTU[id.right]
    n = length(id.right)+1
    seqnum = (max(microbe_y_range)- min( microbe_y_range))/n

    topdat = data.frame(x=  zoom,
                        y = seq(from=min( microbe_y_range),
                                to=max(microbe_y_range),by=seqnum)) %>%
      slice(-1)

    topdat = topdat[-nrow(topdat),]

    topdat$microbe = id.right

    microbe_data_r <- data.frame(
      microbe = topdat$microbe,
      x = topdat$x+ max(heat_coords$x_env  ),
      y = topdat$y
    )

    result_list_r <- list()
    i=1
    for (i in seq_along(microbe_data_r$microbe)) {
      current_microbe <- microbe_data_r$microbe[i]
      r_col <- paste0(current_microbe, "r.BC")
      p_col <- paste0(current_microbe, "p.BC")


      # if (!all(c(r_col, p_col) %in% names(mantel))) next

      result_list_r[[i]] <- mantel %>%
        select(Envs, all_of(c(r_col, p_col))) %>%
        rename(cor = !!r_col, p.value = !!p_col) %>%
        mutate(microbe = current_microbe) %>%
        inner_join(microbe_data_r[i, ], by = "microbe") %>%
        left_join(heat_coords, by = c("Envs" = "Var1")) %>%
        filter(x_env == max(x_env))  # 右侧边缘
      #
      # result_list_r[[i]] <- mantel %>%
      #   select(Envs, all_of(c(r_col, p_col))) %>%
      #   rename(cor = !!r_col, p.value = !!p_col) %>%
      #   mutate(microbe = current_microbe) %>%
      #   inner_join(microbe_data_r[i, ], by = "microbe") %>%
      #   left_join(heat_coords, by = c("Envs" = "Var2")) %>%
      #   distinct(cor.x,.keep_all = T) %>%
      #   mutate( x_env =max(x_env), y_env = 1:length(x_env)  )

    }

    final_df_r <- bind_rows(result_list_r) %>%
      filter(!is.na(x_env), !is.na(y_env)) %>%
      mutate(side = "right")

    final_df_r$x_env = final_df_r$x_env + offset1

    all_dfs$right <- final_df_r
  }



  id.left = matrix.line$left
  if(!is.null(id.left)) {

    #  heatmap_left <- min(heat_coords$x_env) - 1
    microbe_y_range <- range(heat_coords$y_env)
    # tabOTU_l = tabOTU[id.left]
    n = length(id.left)+1
    seqnum = (max(microbe_y_range)- min( microbe_y_range))/n

    topdat = data.frame(x=  -zoom,
                        y = seq(from=min( microbe_y_range),
                                to=max(microbe_y_range),by=seqnum)) %>%
      slice(-1)

    topdat = topdat[-nrow(topdat),]

    topdat$microbe = id.left

    microbe_data_l <- data.frame(
      microbe = topdat$microbe,
      x = topdat$x,
      y = topdat$y
    )

    result_list_l <- list()
    for (i in seq_along(microbe_data_l$microbe)) {
      current_microbe <- microbe_data_l$microbe[i]
      r_col <- paste0(current_microbe, "r.BC")
      p_col <- paste0(current_microbe, "p.BC")

      if (!all(c(r_col, p_col) %in% names(mantel))) next

      result_list_l[[i]] <- mantel %>%
        select(Envs, all_of(c(r_col, p_col))) %>%
        rename(cor = !!r_col, p.value = !!p_col) %>%
        mutate(microbe = current_microbe) %>%
        inner_join(microbe_data_l[i, ], by = "microbe") %>%
        left_join(heat_coords, by = c("Envs" = "Var1")) %>%
        filter(x_env == min(x_env))  # 左侧边缘
      # result_list_l[[i]] <- mantel %>%
      #   select(Envs, all_of(c(r_col, p_col))) %>%
      #   rename(cor = !!r_col, p.value = !!p_col) %>%
      #   mutate(microbe = current_microbe) %>%
      #   inner_join(microbe_data_l[i, ], by = "microbe") %>%
      #   left_join(heat_coords, by = c("Envs" = "Var1")) %>%
      #   distinct(cor.x,.keep_all = T) %>%
      #   mutate( x_env =min(x_env), y_env = 1:length(x_env)  )
      #

    }

    final_df_l <- bind_rows(result_list_l) %>%
      filter(!is.na(x_env), !is.na(y_env)) %>%
      mutate(side = "left")
    final_df_l$x_env = final_df_l$x_env - offset1
    all_dfs$left <- final_df_l


  }

  id.top = matrix.line$up
  if(!is.null(id.top)) {
    #   heatmap_top <- max(heat_coords$y_env) + 1


    microbe_x_range <- range(heat_coords$x_env)

    tabOTU_t = tabOTU[id.top]

    n = length(id.top)+1

    seqnum = (max(microbe_x_range)- min( microbe_x_range))/n

    topdat = data.frame(y=  zoom,
                        x = seq(from=min( microbe_x_range),
                                to=max(microbe_x_range),by=seqnum)) %>%
      slice(-1)

    topdat = topdat[-nrow(topdat),]

    topdat$microbe = id.top

    microbe_data_t <- data.frame(
      microbe = topdat$microbe,
      x = topdat$x,
      y = topdat$y+ max(heat_coords$y_env)
    )


    result_list_t <- list()
    for (i in seq_along(microbe_data_t$microbe)) {
      current_microbe <- microbe_data_t$microbe[i]
      r_col <- paste0(current_microbe, "r.BC")
      p_col <- paste0(current_microbe, "p.BC")

      #  if (!all(c(r_col, p_col) %in% names(mantel))) next

      # result_list_t[[i]] <- mantel %>%
      #   select(Envs, all_of(c(r_col, p_col))) %>%
      #   rename(cor = !!r_col, p.value = !!p_col) %>%
      #   mutate(microbe = current_microbe) %>%
      #   inner_join(microbe_data_t[i, ], by = "microbe") %>%
      #   left_join(heat_coords, by = c("Envs" = "Var1")) %>%
      #   filter(y_env == max(y_env))  # 上侧边缘

      result_list_t[[i]] <- mantel %>%
        select(Envs, all_of(c(r_col, p_col))) %>%
        rename(cor = !!r_col, p.value = !!p_col) %>%
        mutate(microbe = current_microbe) %>%
        inner_join(microbe_data_t[i, ], by = "microbe") %>%
        left_join(heat_coords, by = c("Envs" = "Var1")) %>%
        distinct(cor.x,.keep_all = T) %>%
        mutate(y_env =max(y_env), x_env = 1:length(y_env))

    }

    final_df_t <- bind_rows(result_list_t) %>%
      #  filter(!is.na(x), !is.na(y)) %>%
      mutate(side = "top")
    final_df_t$y_env= final_df_t$y_env + offset1
    all_dfs$top <- final_df_t


  }

  id.bottom = matrix.line$bottom
  if(!is.null(id.bottom)) {


    microbe_x_range <- range(heat_coords$x_env)

    tabOTU_t = tabOTU[id.bottom]

    n = length(id.bottom)+1

    seqnum = (max(microbe_x_range)- min( microbe_x_range))/n

    topdat = data.frame(y=  -zoom,
                        x = seq(from=min( microbe_x_range),
                                to=max(microbe_x_range),by=seqnum)) %>%
      slice(-1)

    topdat = topdat[-nrow(topdat),]

    topdat$microbe = id.bottom

    microbe_data_b <- data.frame(
      microbe = topdat$microbe,
      x = topdat$x,
      y = topdat$y
    )

    result_list_b <- list()
    for (i in seq_along(microbe_data_b$microbe)) {
      current_microbe <- microbe_data_b$microbe[i]
      r_col <- paste0(current_microbe, "r.BC")
      p_col <- paste0(current_microbe, "p.BC")

      if (!all(c(r_col, p_col) %in% names(mantel))) next

      # result_list_b[[i]] <- mantel %>%
      #   select(Envs, all_of(c(r_col, p_col))) %>%
      #   rename(cor = !!r_col, p.value = !!p_col) %>%
      #   mutate(microbe = current_microbe) %>%
      #   inner_join(microbe_data_b[i, ], by = "microbe") %>%
      #   left_join(heat_coords, by = c("Envs" = "Var1")) %>%
      #   filter(y_env == min(y_env))  # 下侧边缘
      #


      result_list_b[[i]] <- mantel %>%
        select(Envs, all_of(c(r_col, p_col))) %>%
        rename(cor = !!r_col, p.value = !!p_col) %>%
        mutate(microbe = current_microbe) %>%
        inner_join(microbe_data_b[i, ], by = "microbe") %>%
        left_join(heat_coords, by = c("Envs" = "Var1")) %>%
        distinct(cor.x,.keep_all = T) %>%
        mutate( y_env =min(y_env), x_env = 1:length(x_env)  )

    }

    final_df_b <- bind_rows(result_list_b) %>%
      #  filter(!is.na(x), !is.na(y)) %>%
      mutate(side = "bottom")
    final_df_b$y_env= final_df_b$y_env - offset1
    all_dfs$bottom <- final_df_b

  }

  final_all <- bind_rows(all_dfs, .id = "direction")

  id1 = (final_all$Envs %>% unique())[1:(floor(length(final_all$Envs %>% unique())/2))]

  id2 = (final_all$Envs %>% unique())[(floor(length(final_all$Envs %>% unique())/2) + 1):
                                        c(length(final_all$Envs %>% unique()))]

  head( final_all)
  group <- c()
  group[ final_all$cor.x > 0] <- "pos"
  group[final_all$cor.x< 0] <- "neg"
  group[final_all$cor.x== 0] <- "nose"
  final_all$groupr = group

  if(mantel_p =="sig"){

    p2 = p +
      # #  ggnewscale::new_scale()+
      # geom_curve(
      #   data = final_all %>% filter(p.value< 0.05),
      #   aes(x = x_env ,
      #       y = y_env ,
      #       xend = x,
      #       yend = y,
      #       color = cor.x,        # 使用相关系数着色
      #       alpha = abs(cor.x),
      #       # curvature = curvature_dir, # 透明度映射相关系数
      #       # angle = angle_adj,
      #       linewidth = abs(cor.x)
      #   ),
      #   curvature = corva,   # 动态曲率
      #   angle = angle ,
      #   lineend = "round",      # 线端形状
      #   show.legend = TRUE
      # ) +
      geom_curve(
        data = final_all %>% filter(p.value< p.thur)%>% filter(Envs %in% id1) %>%filter(side%in% c("right","top")) ,
        aes(x = x_env ,
            y = y_env ,
            xend = x,
            yend = y,
            color = groupr,        # 使用相关系数着色
            # alpha = abs(cor.x),
            # curvature = curvature_dir, # 透明度映射相关系数
            # angle = angle_adj,
            linewidth = abs(cor.x)
        ),
        curvature = corva,   # 动态曲率
        angle =angle ,
        lineend = "round",      # 线端形状
        show.legend = TRUE
      ) +
      geom_curve(
        data = final_all %>% filter(p.value< p.thur)%>% filter(Envs %in% id2) %>%filter(side%in% c("left","bottom")),
        aes(x = x_env ,
            y = y_env ,
            xend = x,
            yend = y,
            color = groupr,        # 使用相关系数着色
            # alpha = abs(cor.x),
            # curvature = curvature_dir, # 透明度映射相关系数
            # angle = angle_adj,
            linewidth = abs(cor.x)
        ),
        curvature = corva,   # 动态曲率
        angle =angle ,
        lineend = "round",      # 线端形状
        show.legend = TRUE
      ) +
      geom_curve(
        data = final_all%>% filter(p.value< p.thur) %>% filter(Envs %in% id1) %>%filter(side%in% c("left","bottom")) ,
        aes(x = x_env ,
            y = y_env ,
            xend = x,
            yend = y,
            color = groupr,        # 使用相关系数着色
            # alpha = abs(cor.x),
            # curvature = curvature_dir, # 透明度映射相关系数
            # angle = angle_adj,
            linewidth = abs(cor.x)
        ),
        curvature = -corva,   # 动态曲率
        angle =angle ,
        lineend = "round",      # 线端形状
        show.legend = TRUE
      ) +
      geom_curve(
        data = final_all %>% filter(p.value< p.thur)%>% filter(Envs %in% id2) %>%filter(side%in% c("right","top")),
        aes(x = x_env ,
            y = y_env ,
            xend = x,
            yend = y,
            color = groupr,        # 使用相关系数着色
            # alpha = abs(cor.x),
            # curvature = curvature_dir, # 透明度映射相关系数
            # angle = angle_adj,
            linewidth = abs(cor.x)
        ),
        curvature = -corva,   # 动态曲率
        angle =angle ,
        lineend = "round",      # 线端形状
        show.legend = TRUE
      ) +
      # # 颜色标度分层设置
      # scale_color_manual(
      #   name = "Microbes",
      #   values = c("#1B9E77", "#D95F02", "#7570B3")  # 手动指定离散颜色
      # ) +
      # scale_color_gradient2(
      #   name = "Correlation",  # 单独为连线设置连续色标
      #   low = "#2E8B57",
      #   mid = "white",
      #   high = "#CD5C5C",
      #   midpoint = 0,
      #   guide = guide_colorbar(
      #     title.position = "top",
      #     barwidth = unit(3, "cm"))
      # ) +
      scale_color_manual(values = c("#91331FFF","#46732EFF")) +
      # scale_fill_manual(values = c("#91331FFF","#46732EFF")) +
      # 透明度和线宽标度
      scale_alpha_continuous(range = c(0.3, 0.8), guide = "none") +
      scale_linewidth_continuous(range = c(0.5, 2)) +
      # 图例布局优化
      guides(
        color = guide_legend(order = 1),              # 微生物颜色图例
        shape = guide_legend(order = 1),              # 形状图例
        linewidth = guide_legend(title = "Cor"), # 线宽图例
        color = guide_colorbar(order = 2)             # 相关图例
      ) +
      # 主题增强
      theme(
        legend.box = "vertical",      # 垂直排列图例
        legend.spacing.y = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm")
      )+
      geom_point(
        data = final_all %>% distinct( x,y, .keep_all = TRUE),
        aes(x = x, y = y),
        shape = 21 , color = "black", fill =  "white",size=3) +
      # 微生物标签）
      ggrepel::geom_text_repel(
        data = final_all %>% distinct( x,y, .keep_all = TRUE),
        aes(x = x , y = y, label = microbe),
        #direction = "y",    # 优先垂直方向调整
        #nudge_x = 0.5,      # 初始水平偏移
        segment.color = NA, # 隐藏连接线
        #min.segment.length = 0,
        # box.padding = 0.2,
        size = 3.5
      )+
      coord_cartesian(
        clip = "off",          # 允许标签溢出绘图区
        expand = FALSE         # 完全禁用扩展
      )
    # ) +
    # theme(
    #   legend.position = c(1 , 0.5),  # 精准定位图例
    #   # legend.box.background = element_rect(fill = "white", color = "black"),
    #   plot.margin = margin(r = 2, unit = "cm",b=1)  # 右侧扩展空间
    #   )

    # # 去掉横纵坐标文字-------
    # # 当四面都有数据时恢复默认
    # if (!is.null(id.top) & !is.null(id.bottom) &
    #     !is.null(id.left) & !is.null(id.right)) {
    #   p2 <- p2 + theme(
    #     axis.text.x = element_text(),
    #     axis.text.y = element_text()
    #   )
    # }
    #
    # # 横坐标处理（上下方向）
    # if (!is.null(id.top) && !is.null(id.bottom)) {
    #   # 情况1：上下都有点 → 隐藏横坐标文字
    #   p2 <- p2 + theme(axis.text.x = element_blank())
    # } else {
    #   # 情况2：至少一侧无点 → 判断具体位置
    #   if (is.null(id.top)) {
    #     # 上边无点 → 坐标轴显示在底部
    #
    #     p2 <- p2 +
    #       scale_x_discrete(position = "top") +
    #       theme(axis.text.x.bottom = element_blank())
    #   }
    #   if (is.null(id.bottom)) {
    #     # 下边无点 → 坐标轴显示在顶部
    #     p2 <- p2 +
    #       scale_x_discrete(position = "bottom") +
    #       theme(axis.text.x.top = element_blank())
    #   }
    # }
    #
    # # 纵坐标处理（左右方向）
    # if (!is.null(id.left) && !is.null(id.right)) {
    #   # 情况1：左右都有点 → 隐藏纵坐标文字
    #   p2 <- p2 + theme(axis.text.y = element_blank())
    # } else {
    #   # 情况2：至少一侧无点 → 判断具体位置
    #   if (is.null(id.left)) {
    #     # 左边无点 → 坐标轴显示在右侧
    #     p2 <- p2 +
    #       scale_y_discrete(position = "left") +
    #       theme(axis.text.y.right = element_blank())
    #   }
    #   if (is.null(id.right)) {
    #     # 右边无点 → 坐标轴显示在左侧
    #     p2 <- p2 +
    #       scale_y_discrete(position = "right") +
    #       theme(axis.text.y.left = element_blank())
    #   }
    # }






  }


  if(mantel_p =="nosig"){
    # 绘图
    p2 = p +
      #  ggnewscale::new_scale()+
      geom_curve(
        data = final_all %>% filter(Envs %in% id1) %>%filter(side%in% c("right","top")) ,
        aes(x = x_env ,
            y = y_env ,
            xend = x,
            yend = y,
            color =groupr,        # 使用相关系数着色
            # alpha = abs(cor.x),
            # curvature = curvature_dir, # 透明度映射相关系数
            # angle = angle_adj,
            linewidth = abs(cor.x)
        ),
        curvature = corva,   # 动态曲率
        angle = angle ,
        lineend = "round",      # 线端形状
        show.legend = TRUE
      ) +
      geom_curve(
        data = final_all %>% filter(Envs %in% id2) %>%filter(side%in% c("left","bottom")),
        aes(x = x_env ,
            y = y_env ,
            xend = x,
            yend = y,
            color = groupr,        # 使用相关系数着色
            # alpha = abs(cor.x),
            # curvature = curvature_dir, # 透明度映射相关系数
            # angle = angle_adj,
            linewidth = abs(cor.x)
        ),
        curvature = corva,   # 动态曲率
        angle =angle ,
        lineend = "round",      # 线端形状
        show.legend = TRUE
      ) +
      geom_curve(
        data = final_all %>% filter(Envs %in% id1) %>%filter(side%in% c("left","bottom")) ,
        aes(x = x_env ,
            y = y_env ,
            xend = x,
            yend = y,
            color = groupr,        # 使用相关系数着色
            # alpha = abs(cor.x),
            # curvature = curvature_dir, # 透明度映射相关系数
            # angle = angle_adj,
            linewidth = abs(cor.x)
        ),
        curvature = -corva,   # 动态曲率
        angle =angle ,
        lineend = "round",      # 线端形状
        show.legend = TRUE
      ) +
      geom_curve(
        data = final_all %>% filter(Envs %in% id2) %>%filter(side%in% c("right","top")),
        aes(x = x_env ,
            y = y_env ,
            xend = x,
            yend = y,
            color = groupr,        # 使用相关系数着色
            # alpha = abs(cor.x),
            # curvature = curvature_dir, # 透明度映射相关系数
            # angle = angle_adj,
            linewidth = abs(cor.x)
        ),
        curvature = -corva,   # 动态曲率
        angle =angle ,
        lineend = "round",      # 线端形状
        show.legend = TRUE
      ) +
      # # 颜色标度分层设置
      # scale_color_manual(
      #   name = "Microbes",
      #   values = c("#1B9E77", "#D95F02", "#7570B3")  # 手动指定离散颜色
      # ) +
      # scale_color_gradient2(
      #   name = "Correlation",  # 单独为连线设置连续色标
      #   low = "#008B45FF",
      #   high = "#EE0000FF",
      #   midpoint = 0,
      #   guide = guide_colorbar(
      #     title.position = "top",
      #     barwidth = unit(3, "cm"))
      # ) +
      scale_color_manual(values = c("#91331FFF","#46732EFF")) +
      # scale_fill_manual(values = c("#91331FFF","#46732EFF")) +
      # 透明度和线宽标度
      # scale_alpha_continuous(range = c(0.3, 0.8), guide = "none") +
      scale_linewidth_continuous(range = c(0.5, 2)) +
      # 图例布局优化
      guides(
        # color = guide_legend(order = 1),              # 微生物颜色图例
        # shape = guide_legend(order = 1),              # 形状图例
        linewidth = guide_legend(title = "Cor"), # 线宽图例
        # color = guide_colorbar(order = 2)             # 相关图例
      ) +
      # 主题增强
      theme(
        legend.box = "vertical",      # 垂直排列图例
        legend.spacing.y = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm")
      )+
      geom_point(
        data = final_all %>% distinct( x,y, .keep_all = TRUE),
        aes(x = x, y = y),
        shape = 21 , color = "black", fill =  "white",size=3) +
      # 微生物标签）
      ggrepel::geom_text_repel(
        data = final_all %>% distinct( x,y, .keep_all = TRUE),
        aes(x = x , y = y, label = microbe),
        #direction = "y",    # 优先垂直方向调整
        #nudge_x = 0.5,      # 初始水平偏移
        segment.color = NA, # 隐藏连接线
        #min.segment.length = 0,
        # box.padding = 0.2,
        size = 3.5
      ) +
      coord_cartesian(
        clip = "off",          # 允许标签溢出绘图区
        expand = FALSE         # 完全禁用扩展
      )
    # theme(
    #   legend.position = c(1 , 0.5),  # 精准定位图例
    #   # legend.box.background = element_rect(fill = "white", color = "black"),
    #   plot.margin = margin(r = 2, unit = "cm",b=1)  # 右侧扩展空间
    #   )

    p2
    # # 去掉横纵坐标文字-------
    # # 当四面都有数据时恢复默认
    # if (!is.null(id.top) & !is.null(id.bottom) &
    #     !is.null(id.left) & !is.null(id.right)) {
    #   p2 <- p2 + theme(
    #     axis.text.x = element_text(),
    #     axis.text.y = element_text()
    #   )
    # }
    #
    # # 横坐标处理（上下方向）
    # if (!is.null(id.top) && !is.null(id.bottom)) {
    #   # 情况1：上下都有点 → 隐藏横坐标文字
    #   p2 <- p2 + theme(axis.text.x = element_blank())
    # } else {
    #   # 情况2：至少一侧无点 → 判断具体位置
    #   if (is.null(id.top)) {
    #     # 上边无点 → 坐标轴显示在底部
    #
    #
    #     p2 <- p2 +
    #       scale_x_discrete(position = "top") +
    #       theme(axis.text.x.bottom = element_blank())
    #   }
    #   if (is.null(id.bottom)) {
    #     # 下边无点 → 坐标轴显示在顶部
    #     p2 <- p2 +
    #       scale_x_discrete(position = "bottom") +
    #       theme(axis.text.x.top = element_blank())
    #   }
    # }
    #
    # # 纵坐标处理（左右方向）
    # if (!is.null(id.left) && !is.null(id.right)) {
    #   # 情况1：左右都有点 → 隐藏纵坐标文字
    #   p2 <- p2 + theme(axis.text.y = element_blank())
    # } else {
    #   # 情况2：至少一侧无点 → 判断具体位置
    #   if (is.null(id.left)) {
    #     # 左边无点 → 坐标轴显示在右侧
    #     p2 <- p2 +
    #       scale_y_discrete(position = "left") +
    #       theme(axis.text.y.right = element_blank())
    #   }
    #   if (is.null(id.right)) {
    #     # 右边无点 → 坐标轴显示在左侧
    #     p2 <- p2 +
    #       scale_y_discrete(position = "right") +
    #       theme(axis.text.y.left = element_blank())
    #   }
    # }
    #

  }

  return(list(p2,final_all))

}




