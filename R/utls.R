#' transform a otu table object from phyloseq object to matrix
#'
#' @param physeq a phyloseq object from R package phyloseq
#' @examples
#' otu = otu_table(ps)
#' @return matrix
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


vegan_otu <-  function(physeq){
  OTU <-  otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <-  t(OTU)
  }
  return(as(OTU,"matrix"))
}


#' transform a tax object from phyloseq object to matrix
#'
#' @param physeq a phyloseq object from R package phyloseq
#' @examples
#' tax = vegan_tax(ps)
#' @return matrix
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export

vegan_tax <-  function(physeq){
  tax <-  tax_table(physeq)

  return(as(tax,"matrix"))
}


#' Quickly import microbiome data and bulding phyloseq object
#'
#' @param otu dataframe,otu table
#' @param tax data.frame, tax table
#' @param map matrix or dataframe, including sampleID and groupID;
#' @param tree Evolutionary tree
#' @param ps a phyloseq object
#' @param group column name for groupID.
#' @examples
#' data(otu)
#' data(tax)
#' data(map)
#' inputMicro(otu,tax,map,tree,group  = group)
#' @return a phyloseq object
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export





inputMicro = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group"){

  if (is.null(otu)&is.null(tax)&is.null(map)) {
    ps = ps
    map = as.data.frame(sample_data(ps))
    map = map[, group]
    colnames(map) = "Group"
    map$Group = as.factor(map$Group)
    sample_data(ps) = map
    map = NULL
  }

  if (is.null(ps) ) {

    if (!is.null(otu)) {
      head(otu)
      otu = as.matrix(otu)
      str(otu)

      ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE))

    }

    if (!is.null(tax) ) {
      head(tax)
      tax = as.matrix(tax)
      # taxa_names(tax)
      x = tax_table(tax)
      ps = merge_phyloseq(ps,x)
      ps
    }


    if (!is.null(map) ){

      map = map[group]

      map[,group] = as.factor(map[,group] )
      map$Group
      z  = sample_data(map)
      ps = merge_phyloseq(ps,z)
      ps
    }
    if (!is.null(tree) ) {
      # #导入进化树
      h = phy_tree(tree)
      ps = merge_phyloseq(ps,h)
      ps
    }


  }
  return(ps)

}



