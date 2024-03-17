# 
# #  模块比对，定义相似模块#--------
# res1 = module.compare.m(
#   ps = NULL,
#   corg = cortab,
#   zipi = FALSE,
#   zoom = 0.2,
#   padj = F,
#   n = 3)
# 
# #不同分组使用一个圆圈展示，圆圈内一个点代表一个模块，相连接的模块代表了相似的模块。p1 = res1[[1]]
# # p1
# #--提取模块的OTU，分组等的对应信息
# dat1 = res1[[2]]
# head(dat1)
# colnames(dat1)[2] = "group.module"
# #模块相似度结果表格
# dat2 = res1[[3]]
# head(dat2)
# head(node)
# 
# 
# node4 = data.module.Facet.line (
#   dat1 = dat1,#模块的OTU，分组等的对应信息
#   dat2 = dat2,#模块相似度结果表格
#   node = node, #包含节点坐标的节点
#   n = 12# 展示丰度最高的n个模块
# )


#  使用node节点表格和module.compare.m函数的输出内容
data.module.Facet.line = function(
    node  = node, #包含节点坐标的节点
    dat1,#模块的OTU，分组等的对应信息
    dat2,#模块相似度结果表格
    n = 12# 展示丰度最高的n个模块
){
  #  拆分模块
  dat2$module2 = row.names(dat2) %>% strsplit("[-]") %>%
    sapply(`[`, 2)
  dat2$module1 = row.names(dat2) %>% strsplit("[-]") %>%
    sapply(`[`, 1)
  dat2$m1 = dat2$module1 %>% strsplit("model") %>%
    sapply(`[`, 1)
  dat2$m2 = dat2$module2 %>% strsplit("model") %>%
    sapply(`[`, 1)
  dat2$cross = paste(dat2$m1,dat2$m2,sep = "_Vs_")
  # head(dat2)
  dat2 = dat2 %>% filter(module1 != "none")
  
  
  
  # 相似模块合并#------
  id = dat2$module1 %>% unique()
  for (i in 1:length(id)) {
    tem1 = dat2 %>% 
      filter(module1 == id[i]) %>% 
      filter(p_adj < 0.05) %>% .$module2
    tem2 = c(id[i],tem1)
    tem3 = data.frame(ID = tem2,class = paste0("sameM.",i))
    
    if (i == 1) {
      tem4 = tem3
    } else{
      tem4 = rbind(tem4,tem3)
    }
  }
  
  tem4$ID %>% 
    unique() %>%
    length()
  
  id = tem4$ID %>% table() %>% as.data.frame() %>% 
    arrange(desc(Freq)) %>%
    filter(Freq > 1) %>%
    .$`.` %>% as.character()
  
  
  for (i in 1:length(id)) {
    tem5 = tem4 %>% filter(ID == id[i]) %>%.$class
    print( tem5)
    
    sid = tem5[str_detect(tem5,"MNEW")]
    if (length(sid)> 0) {
      for (ii in 1:length(tem5)) {
        tem4$class = gsub(paste0(tem5[ii],"$"),sid,tem4$class,perl =TRUE)
      }
    } else{
      for (ii in 1:length(tem5)) {
        tem4$class = gsub(paste0(tem5[ii],"$"),paste0("MNEW",i),tem4$class,perl =TRUE)
      }
    }
    tem4 = tem4 %>% distinct(ID,class,.keep_all = TRUE)
  }
  
  tem4$class %>% unique()
  
  
  #  将合并结果加入到节点中
  node3 = node %>% left_join(dat1,by = c("ID","Group"))
  head(node3)
  node4 = node3 %>% left_join(tem4,by = c("group.module" = "ID"))
  head(node4)
  node4  = add.id.facet(node4,"Group")
  head(node4)
  tt.1 = node4$group.module %>% table() %>% as.data.frame() %>% arrange(desc(Freq)) %>%.$`.` %>%
    head(n) %>% as.character()
  # [1:3]
  node4$label = factor(node4$label,levels = as.character(unique( node4$label)))
  node4$Group = factor(node4$Group,levels = node4$Group %>% unique())
  head(node4)
  return(node4)
}
