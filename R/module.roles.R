#'Module roles
#'
#' This function assigns roles to features identified in sub communities using
#' two metrics that is: within-module degree which measures how well a particular feature is connected to others in the same subcommunity (module.)
#' The second metric is among-module connectivity which measures how a feature is linked
#' to other modules in the network. Features are classified as ultra peripherals,
#' peripherals, provincial, connectors, kinless, module hubs, or non hubs.
#'
#' @param comm_graph : Graph object returned by `co_occurence_network` function.
#'
#' @examples
#' taxa.roles <- module.roles(g)
#' p <- plot_roles(taxa.roles)
#' print(p)
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#' @references Guimera, Roger, and Luis A Nunes Amaral. 2005. “Functional Cartography of Complex Metabolic Networks.”
#' Nature 433 (7028). NIH Public Access: 895.
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},  Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export module.roles
#' @export plot_roles


module.roles <- function(comm_graph){

  td <- network_degree(comm_graph)

  wmd <- within_module_degree(comm_graph)

  z <- zscore(wmd)

  amd <- among_module_connectivity(comm_graph)

  pc <- participation_coeffiecient(amd, td)

  zp <- data.frame(z,pc)

  nod.roles <- assign_module_roles(zp)

  return(nod.roles)

}

# find total degree for each of the features in the graph

network_degree <- function(comm_graph){

  ki_total <-NULL

  net_degree <- igraph::degree(comm_graph)

  for(i in 1:length(V(comm_graph))){

    ki <- net_degree[i]

    tmp <- data.frame(taxa=names(ki), total_links=ki)

    if(is.null(ki_total)){ki_total<-tmp} else{ki_total <- rbind(ki_total, tmp)}

  }

  return(ki_total)

}


#compute within-module degree for each of the features

within_module_degree <- function(comm_graph){

  mods <- igraph::get.vertex.attribute(comm_graph, "module")

  vs <- as.list(V(comm_graph))

  modvs <- data.frame("taxon"= names(vs), "mod"=mods)

  sg1 <- decompose.graph(comm_graph,mode="strong")

  df <- data.frame()

  for(mod in unique(modvs$mod)){

    mod_nodes <- subset(modvs$taxon,modvs$mod==mod)

    neighverts <- unique(unlist(sapply(sg1,FUN=function(s){if(any(V(s)$name %in% mod_nodes)) V(s)$name else NULL})))

    g3 <- induced.subgraph(graph=comm_graph,vids=neighverts)

    mod_degree <- igraph::degree(g3)

    for(i in mod_nodes){

      ki <- mod_degree[which(names(mod_degree)==i)]

      tmp <- data.frame(module=mod, taxa=names(ki), mod_links=ki)

      df <- rbind(df,tmp)

    }

  }

  return(df)

}
#calculate the degree (links) of each node to nodes in other modules.

among_module_connectivity <- function(comm_graph){

  mods <- igraph::get.vertex.attribute(comm_graph, "module")

  vs <- as.list(V(comm_graph))

  modvs <- data.frame("taxa"= names(vs), "mod"=mods)

  df <- data.frame()

  for(i in modvs$taxa){

    for(j in modvs$taxa){

      if(are_adjacent(graph=comm_graph, v1=i , v2=j)){

        mod <- subset(modvs$mod, modvs$taxa==j)

        tmp <- data.frame(taxa=i, taxa2=j, deg=1, mod_links=mod)

        df <- rbind(df, tmp)

      }

    }

  }

  out <- aggregate(list(mod_links=df$deg), by=list(taxa=df$taxa, module=df$mod), FUN=sum)

  return(out)

}

#compute within-module degree z-score which
#measures how well-connected a node is to other nodes in the module.

zscore <- function(mod.degree){

  ksi_bar <- aggregate(mod_links ~ module, data=mod.degree, FUN = mean)

  ksi_sigma <- aggregate(mod_links ~ module, data=mod.degree, FUN = sd)

  z <- NULL

  for(i in 1:dim(mod.degree)[1]){

    mod_mean <- ksi_bar$mod_links[which(ksi_bar$module == mod.degree$module[i])]

    mod_sig <- ksi_sigma$mod_links[which(ksi_bar$module == mod.degree$module[i])]

    z[i] <- (mod.degree$mod_links[i] - mod_mean)/mod_sig

  }

  z <- data.frame(row.names=rownames(mod.degree), z, module=mod.degree$module)

  return(z)

}


#The participation coefficient of a node measures how well a  node is distributed
# in the entire network. It is close to 1 if its links are uniformly
#distributed among all the modules and 0 if all its links are within its own module.

participation_coeffiecient <- function(mod.degree, total.degree){

  p <- NULL

  for(i in total.degree$taxa){

    ki <- subset(total.degree$total_links, total.degree$taxa==i)

    taxa.mod.degree <- subset(mod.degree$mod_links, mod.degree$taxa==i)

    p[i] <- 1 - (sum((taxa.mod.degree)**2)/ki**2)

  }

  p <- as.data.frame(p)

  return(p)

}


assign_module_roles <- function(zp){

  zp <- na.omit(zp)

  zp$roles <- rep(0, dim(zp)[1])

  outdf <- NULL

  for(i in 1:dim(zp)[1]){

    df <- zp[i, ]

    if(df$z < 2.5){ #non hubs

      if(df$p < 0.05){

        df$roles <- "ultra peripheral"

      }
      else if(df$p < 0.620){

        df$roles <- "peripheral"

      }
      else if(df$p < 0.80){

        df$roles <- "non hub connector"

      }
      else{

        df$roles <- "non hub kinless"

      }

    }
    else { # module hubs

      if(df$p < 0.3){

        df$roles <- "provincial hub"

      }
      else if(df$p < 0.75){

        df$roles <- "connector hub"

      }
      else {

        df$roles <- "kinless hub"

      }

    }

    if(is.null(outdf)){outdf <- df}else{outdf <- rbind(outdf, df)}

  }

  return(outdf)

}

plot_roles <- function(node.roles, roles.colors=NULL){

  x1<- c(0, 0.05, 0.62, 0.8, 0, 0.30, 0.75)
  x2<- c(0.05, 0.62, 0.80, 1,  0.30, 0.75, 1)
  y1<- c(-Inf,-Inf, -Inf, -Inf,  2.5, 2.5, 2.5)
  y2 <- c(2.5,2.5, 2.5, 2.5, Inf, Inf, Inf)

  lab <- c("ultra peripheral","peripheral" ,"non-hub connector","non-hub kinless","provincial"," hub connector","hub kinless")

  if(is.null(roles.colors)){roles.colors <- c("#E6E6FA", "#DCDCDC", "#F5FFFA", "#FAEBD7", "#EEE8AA", "#E0FFFF", "#F5F5DC")}

  p <- ggplot() + geom_rect(data=NULL, mapping=aes(xmin=x1, xmax=x2, ymin=y1,ymax=y2, fill=lab))

  p <- p + guides(fill=guide_legend(title="Topological roles"))

  p  <- p + scale_fill_manual(values = roles.colors)

  p <- p + geom_point(data=node.roles, aes(x=p, y=z,color=module)) + theme_bw()

  p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("Participation Coefficient")+ylab(" Within-module connectivity z-score")

  return(p)
}

