merge16S_ITS <- function(ps16s = ps16s,
                         psITS = psITS,
                         N16s = 100,
                         NITS = 100,
                         scale = TRUE,
                         onlygroup = FALSE,
                         dat1.lab = "bac",
                         dat2.lab = "fun"
) {


  if (scale == TRUE) {
    if (!is.null(ps16s)) {
      ps16s  = phyloseq::transform_sample_counts(ps16s, function(x) x / sum(x) )
    }
    if (!is.null(psITS)) {
      psITS  = phyloseq::transform_sample_counts(psITS, function(x) x / sum(x) )
    }

  }
  if (!is.null(ps16s)) {
    # ps_16s = phyloseq::filter_taxa(ps16s, function(x) mean(x) > N16s, TRUE)#select OTUs according to  relative abundance
    ps_16s  =  filter_OTU_ps(ps = ps16s,Top = N16s)
    ###
    otu_table_16s = as.data.frame(t(vegan_otu(ps_16s)))
    row.names(otu_table_16s) = paste(dat1.lab,row.names(otu_table_16s),sep = "_")
    ## change the OTU name of bac and fungi OTU table
    tax_table_16s = as.data.frame(vegan_tax(ps_16s))
    row.names(tax_table_16s) = paste(dat1.lab,row.names(tax_table_16s),sep = "_")
    #-- add a col marked the bac and fungi
    tax_table_16s$filed = rep(dat1.lab,length(row.names(tax_table_16s)))
  }
  if (!is.null(psITS)) {
    ps_ITS = phyloseq::filter_taxa(psITS, function(x) mean(x) > NITS , TRUE)#select OTUs according to  relative abundance
    ps_ITS = filter_OTU_ps(ps = psITS,Top = NIS)
    otu_table_ITS = as.data.frame(t(vegan_otu(ps_ITS)))
    row.names(otu_table_ITS) = paste(dat2.lab,row.names(otu_table_ITS ),sep = "_")
    tax_table_ITS = as.data.frame(vegan_tax(ps_ITS))
    row.names(tax_table_ITS) = paste(dat2.lab,row.names(tax_table_ITS),sep = "_")
    tax_table_ITS$filed = rep(dat2.lab,length(row.names(tax_table_ITS)))

  }


  if (!is.null(psITS) & !is.null(ps16s) ) {
    ## merge OTU table of bac and fungi
    otu_table = rbind(otu_table_16s,otu_table_ITS)

    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_16s,tax_table_ITS)
      dim(otu_table)
    } else if(onlygroup == TRUE){
      tax_table = data.frame(filed = c(tax_table_16s$filed,tax_table_ITS$filed),row.names = row.names(otu_table),id = row.names(otu_table))
    }
    #on of map table as final map table
    mapping = as.data.frame(sample_data(ps_16s))
    head(mapping)
    # mapping$Group4 = "all_sample"
    # mapping$Group4 = as.factor(mapping$Group4)
    ##merge all abject of phyloseq
    pallps <- phyloseq(otu_table(as.matrix(otu_table),taxa_are_rows = T),
                       sample_data(mapping),
                       tax_table(as.matrix(tax_table)))


  } else if(is.null(psITS) & !is.null(ps16s) ) {
    otu_table = rbind(otu_table_16s)

    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_16s)
      dim(otu_table)
    } else if(onlygroup == TRUE){
      tax_table = data.frame(filed = c(tax_table_16s$filed,row.names = row.names(otu_table),id = row.names(otu_table)))
    }
    #on of map table as final map table
    mapping = as.data.frame(sample_data(ps_16s))
    head(mapping)
    # mapping$Group4 = "all_sample"
    # mapping$Group4 = as.factor(mapping$Group4)
    ##merge all abject of phyloseq
    pallps <- phyloseq(otu_table(as.matrix(otu_table),taxa_are_rows = T),
                       sample_data(mapping),
                       tax_table(as.matrix(tax_table)))


  } else if (!is.null(psITS) & is.null(ps16s)){
    otu_table = rbind(otu_table_ITS)

    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_ITS)
      dim(otu_table)
    } else if(onlygroup == TRUE){
      tax_table = data.frame(filed = c(tax_table_ITS$filed),row.names = row.names(otu_table),id = row.names(otu_table))
    }
    #on of map table as final map table
    mapping = as.data.frame(sample_data(psITS))
    head(mapping)
    # mapping$Group4 = "all_sample"
    # mapping$Group4 = as.factor(mapping$Group4)
    ##merge all abject of phyloseq
    pallps <- phyloseq(otu_table(as.matrix(otu_table),taxa_are_rows = T),
                       sample_data(mapping),
                       tax_table(as.matrix(tax_table)))

  }
  return(pallps)
}
