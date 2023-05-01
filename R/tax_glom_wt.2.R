

# 如果发现OTU名字中有特殊字符，无法作为行名，则自动转化为一个list，包含OTU表格和注释文件。
tax_glom_wt.2 <- function(ps = ps,
                        ranks = "Phylum") {
  
  
  if (  is.numeric(ranks)) {
    ranks <- phyloseq::rank.names(ps)[ranks]
  }
  
  
  otu <- as.data.frame(t(vegan_otu(ps)))
  tax <- as.data.frame(vegan_tax(ps))
  
  # building group
  tax[[ranks]][is.na(tax[[ranks]])] = "Unknown"
  tax[[ranks]][tax[[ranks]] == ""] = "Unknown"
  tax[[ranks]][tax[[ranks]] == "NA"] = "Unknown"
  split <- split(otu,tax[[ranks]])
  #calculate sum by group
  apply <- lapply(split,function(x)colSums(x[]))
  # result chack
  otucon <- do.call(rbind,apply)
  
  taxcon <- tax[1:match(ranks,colnames(tax))]
  taxcon <- taxcon[!duplicated(tax[[ranks]]),]
  
  if (is.vector(taxcon)) {
    taxcon = data.frame(row.names = taxcon,ranks = taxcon)
    colnames(taxcon) = ranks
  }
  
  #-tax name with NA wound be repeated with unknown
  taxcon[[ranks]][is.na(taxcon[[ranks]])] = "unknown"
  row.names(taxcon) <- taxcon[[ranks]]
  
  

  
  fit<-try(
    pscon <- phyloseq::phyloseq(
      phyloseq::otu_table( as.matrix(otucon),taxa_are_rows = TRUE),
      phyloseq::tax_table(as.matrix(taxcon)),
      phyloseq::sample_data(ps)
    )
  )
  if("try-error" %in% class(fit))
  {
    pscon = list(otucon = otucon,
                 taxcon = taxcon)
    
  }else{
   
  }
  
  
  
  return(pscon)
}
