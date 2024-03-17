Facet.line.plot = function(
    g = p,
    node4 = node4,# 节点表格
    id.facet = id.facet, # 分面坐标排序
    class= class,# 模块信息
    group = group# 分面的分组信息
){
  node4$group = as.factor(node4$group)
  id = node4$class %>% unique()
  id = id[!is.na(id)]
  for (j in 1:length(id)) {
    # node4 %>% filter(class == id[j]) %>% distinct(class,group,.keep_all = TRUE)
    iid = node4 %>% filter(class == !!syms(id[j])) %>% 
      distinct(class,group,.keep_all = TRUE) %>% .$id.facet
    iid2 = iid %>% strsplit( "[_]") %>% sapply(`[`, 2)  %>% 
      as.numeric()
    # iid2 = c(1,1,2)
    iid3 = iid %>% strsplit( "[_]") %>% sapply(`[`, 1) 
    # iid3 = c("Group3","Group1")
    iid3 = levels(node4$group) [levels(node4$group) %in% iid3]
    iid4 = (1:length(levels(node4$group)))[levels(node4$group) %in% iid3]
    tab = data.frame(Group =  iid3,facet = iid4)
    tab2 = tab$facet[tab$Group%in%c(node4 %>% filter(class == !!syms(id[j])) %>% 
                                      distinct(class,group,.keep_all = TRUE) %>% .$group)]
    
    
    tem = data.frame(ID = c(0,iid2[1:(length(iid2)-1)]),
                     ID2 = iid2,
                     facet =c(0, tab2[1:(length(tab2)-1)]),
                     tab2) 
    tem2 = tem [-1,]
    tem2
    for (i in 1:dim(tem2)[1]) {
      g <- line.across.facets.network(g, 
                                      from=tem2[i,3], to=tem2[i,4], 
                                      from_point_id=tem2[i,1],
                                      to_point_id=tem2[i,2])
      
    }
  }
  return(g)
}
