
#' Hyperbolic distance between points
#'
#' Computes the hyperbolic distance between a point at polar coordinates (\code{ZI$r}, \code{ZI$theta}) and \code{m} points at polar coordinates (\code{ZJ$r}, \code{ZJ$theta}).
#' 
#' @param ZI single-element data frame; A node in hyperbolic space placed at polar coordinates (\code{ZI$r}, \code{ZI$theta}).
#' @param ZJ data frame; m entries representing m nodes in hyperbolic space placed at polar coordinates (\code{ZJ$r}, \code{ZJ$theta}).
#' 
#' @return An \code{m}-element vector with the hyperbolic distance between node \code{ZI} and nodes \code{ZJ}.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Krioukov, D. et al. (2010) Hyperbolic geometry of complex networks. \emph{Physical Review E} 82(3).
#' 
#' @examples
#' # Generate an artificial network with the PS model
#' net <- ps_model(500, 6, 2.5, 0)
#' # If nodes data frame contains m nodes, use the following 
#' # to compute the distance between node 1 and the rest (including itself):
#' nodes <- net$polar
#' d <- hyperbolic_dist(nodes[1,], nodes)
#' # To compute the pairwise distances between m nodes use:
#' D <- sapply(seq(nrow(nodes)), function(i) hyperbolic_dist(nodes[i, ], nodes))
#' 
#' @export
#' 
hyperbolic_dist <- function(ZI, ZJ){
  
  delta <- pi - abs(pi - abs(ZI$theta - ZJ$theta)) # Angular separation between points
  
  d <- cosh(ZI$r)*cosh(ZJ$r) - sinh(ZI$r)*sinh(ZJ$r)*cos(delta)
  
  # Due to precision problems, numbers d < 1 should be set to 1 to get the right final hyperbolic distance of 0
  d[d < 1] <- 1
  
  d <- acosh(d)
  
  #In addition, if ZI == ZJ, d should be 0
  d[ZI$r == ZJ$r & ZI$theta == ZJ$theta] <- 0
  
  return(d)
}

#' Plot a network in hyperbolic space
#' 
#' Plots a network in hyperbolic space, given node polar coordinates and colours.
#' 
#' @param network igraph; The network to be plotted.
#' @param nodes data frame; Polar coordinates (r, theta) of all nodes in the network.
#' @param node.colour vector; An optional parameter that allows to colour nodes according to a colour vector or a single colour.
#' 
#' @return A ggplot of the given network in hyperbolic space.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Generate an artificial network with the PS model
#' net <- ps_model(500, 6, 2.5, 0)
#' # Plot the network and colour nodes according to their radial coordinates
#' plot_hyperbolic_net(net$network, net$polar, net$polar$theta)
#' 
#' @export
#' @import igraph
#' @import ggplot2
#' 
plot_hyperbolic_net <- function(network, nodes, node.colour = 1){
  
  # Change node coordinates from polar to cartesian
  nodes.cart <- data.frame(x = nodes$r * cos(nodes$theta), y = nodes$r * sin(nodes$theta), theta = nodes$theta, node.colour = factor(node.colour))
  edges <- as_edgelist(network, names = F)
  edge.df <- data.frame(x1 = nodes.cart$x[edges[, 1]], y1 = nodes.cart$y[edges[, 1]], x2 = nodes.cart$x[edges[, 2]], y2 = nodes.cart$y[edges[, 2]])
  
  return(ggplot(edge.df, aes_(~x1, ~y1)) + geom_segment(aes_(xend = ~x2, yend = ~y2), size = 0.1, colour="grey", alpha = 0.9) + 
    geom_point(data = nodes.cart, aes_(~x, ~y, colour = ~node.colour), size = 3) + 
    scale_colour_gradientn(colours = grDevices::rainbow(4)[4:1]) + theme_bw() + 
    theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
          axis.text.x = element_blank(), axis.text.y = element_blank(),  
          axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) + 
    labs(title = "Network in the hyperbolic disc", colour = expression(theta)))
}

#' Plot the nodes of a network in hyperbolic space
#' 
#' Plots nodes of a network in hyperbolic space, given node polar coordinates and colours.
#' 
#' @param nodes data frame;  Polar coordinates (r, theta) of all nodes in the network.
#' @param node.colour vector; An optional parameter that allows to colour nodes according to a colour vector or a single colour.
#' 
#' @return A ggplot of the given network nodes in hyperbolic space.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Generate an artificial network with the PS model
#' net <- ps_model(500, 6, 2.5, 0)
#' # Plot the network nodes and colour them all in blue
#' plot_hyperbolic_nodes(net$polar, "blue")
#' 
#' @export
#' @import ggplot2
#' 
plot_hyperbolic_nodes <- function(nodes, node.colour){
  # Change node coordinates from polar to cartesian
  nodes.cart <- data.frame(x = nodes$r * cos(nodes$theta), y = nodes$r * sin(nodes$theta), node.colour = factor(node.colour))
  
  return(ggplot(nodes.cart, aes_(~x, ~y, colour = ~node.colour), size = 3) + geom_point() + scale_colour_manual(values = unique(node.colour)) + 
           theme_bw() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
                              axis.text.x = element_blank(), axis.text.y = element_blank(),  
                              axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) + 
    labs(title = "Network in the hyperbolic disc (nodes only)") + guides(colour = F))
}

#' Generate a log10-linearly spaced sequence
#' 
#' Generates a sequence of \code{n} logarithmically spaced points between \code{s} and \code{e}, inclusive.
#' 
#' @param s numeric; The sequence starting point.
#' @param e numeric; The sequence ending point.
#' @param n integer; The number of points to be generated.
#' 
#' @return A numeric vector with \code{n} logarithmically spaced points between \code{s} and \code{e}.
#' 
log_seq <- function(s, e, n){
  s <- log10(s)
  e <- log10(e)
  return(exp(log(10) * seq(s, e, length.out = n)))
}

#' Plots the degree distribution of a network in log-log scale.
#' 
#' @param network igraph; The network whose degree distribution is to be plotted.
#' @param bins integer; The number of bins.
#' 
#' @return A ggplot of the network's degree distribution.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples 
#' # Generate an artificial network with the PS model
#' net <- ps_model(500, 6, 2.5, 0)
#' # Plot its degree distribution
#' plot_degree_distr(net$network, bins = 50)
#' 
#' @export
#' @import igraph
#' @import ggplot2
#' @importFrom scales trans_breaks trans_format math_format
#' 
plot_degree_distr <- function(network, bins = 100){
  # Calculate degrees for each node
  d <- degree(network, mode = "all")
  
  # Bin the degrees using logarithmic spacing
  breaks = log_seq(min(d), max(d), bins + 1)
  cuts <- cut(d, breaks = breaks, labels = breaks[-length(breaks)], include.lowest = T)
  stats <- data.frame(deg = as.numeric(levels(cuts)), 
                       prob = as.numeric(table(cuts))/sum(as.numeric(table(cuts))), 
                       stringsAsFactors = F)
  stats <- stats[stats$prob > 0, ]
  
  return(ggplot(stats, aes_(~deg, ~prob)) + geom_point() + 
           scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                         labels = scales::trans_format("log10", scales::math_format())) + 
           scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                         labels = scales::trans_format("log10", scales::math_format())) + 
           annotation_logticks() + labs(x = "Node degree", y = "Probability") + theme_bw() +
           theme(panel.grid.minor = element_blank()))
}

#' Check if a given network is strongly clustered
#' 
#' Compares the clustering of a given network with that of \code{epochs} ER graphs with the same number of nodes and edges.
#' 
#' @param network igraph; The network of interest.
#' @param epochs integer; The number of ER graphs to consider in the comparison.
#' @param label character; A label for the given network in the resulting plot.
#' 
#' @return A ggplot object comparing the clustering coefficient of the given network with the average of \code{epochs} ER graphs.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples 
#' # Generate an artificial network with the PS model
#' net <- ps_model(500, 6, 2.5, 0)
#' # Compare its clustering coefficient with ER networks 
#' # with the same number of nodes and edges
#' compare_clustering_to_er(net$network, 100, "PS network")
#' 
#' @export
#' @import igraph
#' @import ggplot2
#' 
compare_clustering_to_er <- function(network, epochs = 100, label = "Given network"){
  er.clust <- numeric(length = epochs)
  N <- vcount(network)
  E <- vcount(network)
  for(i in 1:epochs){
    g <- sample_gnm(N, E, directed = F, loops = F)
    er.clust[i] <- transitivity(g, type = "average")
  }
  df <- data.frame(net = factor(c(label, "ER graphs"), levels = c(label, "ER graphs"), ordered = T), 
                   clust = c(transitivity(network, type = "average"), mean(er.clust)),
                   err = c(0, stats::sd(er.clust)))
  dodge <- position_dodge(width = 0.9)
  return(ggplot(df, aes_(~net, ~clust)) + geom_bar(position = dodge, stat = "identity") + 
           geom_errorbar(aes(ymin = df$clust - df$err, ymax = df$clust + df$err), position = dodge, width = 0.25) + 
           labs(x = "", y = "Clustering coefficient") + theme_bw())
}

#' Analyse the navigability of a complex network
#' 
#' Given a network and source and destination nodes, greedy route packets in hyperbolic space and record successful deliveries. If the packet
#' visits any of the nodes in \code{faulty}, it is dropped and the delivery flagged as unsuccessful (unless the faulty node is a target).
#' 
#' @param net igraph; A complex network with N nodes.
#' @param polar data frame; Polar coordinates (r, theta) of all nodes in the network.
#' @param source vector; A vector with one or more node indices representing sources.
#' @param target vector; A vector, the same size as \code{source}, with one or more node indices representing targets.
#' @param faulty vector; A vector with one or more node indices representing faulty system components (default = c()).
#' 
#' @return A vector, the same size as \code{source} or \code{target}, with the hop stretch required to deliver the packet 
#' from each source to each target, or 0 if the packet was dropped.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Boguna, M. and Krioukov, D. (2009) Navigating Ultrasmall Worlds in Ultrashort Time. \emph{Physical Review Letters} 102(5).
#' @references Krioukov, D. et al. (2010) Hyperbolic geometry of complex networks. \emph{Physical Review E} 82(3).
#' 
#' @examples 
#' # Generate an artificial network with the PS model
#' N <- 500
#' net <- ps_model(N, 6, 2.5, 0)
#' 
#' # Form vectors of random non-redundant source-target pairs
#' st <- 1000
#' # We subtract 1, because the formulae to go from linear upper 
#' # diagonal indexing to (i,j) are zero-based
#' k <- sample(N*(N-1)/2, st) - 1
#' sources <- (N - 2 - floor(sqrt(-8*k + 4*N*(N-1)-7)/2.0 - 0.5))
#' targets <- (k + sources + 1 - N*(N-1)/2 + (N-sources)*((N-sources)-1)/2)
#' 
#' # Back to 1-based indexing
#' sources <- sources + 1
#' targets <- targets + 1
#' 
#' # Analyse the network's navigability
#' hop.stretch <- greedy_route_packets(net$network, net$polar, sources, targets)
#' 
#' # Compute the fraction of succesfully delivered packets
#' sum(hop.stretch > 0)/st
#' 
#' @export
#' @import igraph
#' 
greedy_route_packets <- function(net, polar, source, target, faulty = c()){
  hops <- vector("numeric", length = length(source)) #If routing is successful, hops > 0
  sp <- vector("numeric", length = length(source))
  
  for(i in 1:length(source)){
    status <- "in progress"
    sp[i] <- distances(net, v = source[i], to = target[i], weights = NA)
    src <- source[i]
    previous <- 0
    while(status == "in progress"){
      if(src %in% faulty){
        hops[i] <- 0
        status <- "dropped"
      }
      else{
        neighs <- neighbors(net, v = src)
        closest.to.target <- neighs[which.min(hyperbolic_dist(polar[target[i], ], polar[neighs, ]))]
        if(closest.to.target == target[i]){
          hops[i] <- hops[i] + 1
          status <- "done"
        }else if(closest.to.target == previous){
          hops[i] <- 0
          status <- "dropped"
        }else{
          previous <- src
          src <- closest.to.target
          hops[i] <- hops[i] + 1
        }
      }
    }
  }
  
  return(hops/sp)
}

#' Obtain the greedy paths of a greedy routing process
#' 
#' Given a network and source and destination nodes, greedy route packets in hyperbolic space and record the path followed until the target is reached. 
#' If the packet visits any of the nodes in the vector \code{faulty}, the second closest neighbour to the target is used instead, unless the faulty node is 
#' a target or the source itself.
#' If the packet cannot reach its destination, a \code{-1} is added to the path.
#' 
#' @param net igraph; A complex network with N nodes.
#' @param polar data frame; Polar coordinates (r, theta) of all nodes in the network.
#' @param source vector; A vector with one or more node indices representing sources.
#' @param target vector; A vector, the same size as \code{source}, with one or more node indices representing targets.
#' @param faulty vector; A vector with one or more node indices representing faulty system components (default = c()).
#' 
#' @return A list, the same size of \code{source} or \code{target}, with the paths followed to deliver the packet from each source to each target.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Boguna, M. and Krioukov, D. (2009) Navigating Ultrasmall Worlds in Ultrashort Time. \emph{Physical Review Letters} 102(5).
#' @references Krioukov, D. et al. (2010) Hyperbolic geometry of complex networks. \emph{Physical Review E} 82(3).
#' 
#' @examples 
#' # Generate an artificial network with the PS model
#' N <- 500
#' net <- ps_model(N, 6, 2.5, 0)
#' 
#' # Form vectors of random non-redundant source-target pairs
#' st <- 1000
#' # We subtract 1, because the formulae to go from linear upper 
#' # diagonal indexing to (i,j) are zero-based
#' k <- sample(N*(N-1)/2, st) - 1
#' sources <- (N - 2 - floor(sqrt(-8*k + 4*N*(N-1)-7)/2.0 - 0.5))
#' targets <- (k + sources + 1 - N*(N-1)/2 + (N-sources)*((N-sources)-1)/2)
#' 
#' # Back to 1-based indexing
#' sources <- sources + 1
#' targets <- targets + 1
#' 
#' # Analyse the network's navigability
#' paths <- get_greedy_routing_paths(net$network, net$polar, sources, targets)
#' 
#' @export
#' @import igraph
#' 
get_greedy_routing_paths <- function(net, polar, source, target, faulty = c()){
  paths <- vector("list", length = length(source))
  
  for(i in 1:length(source)){
    status <- "in progress"
    src <- source[i]
    previous <- 0
    p <- src
    if(src %in% faulty){
      status <- "dropped"
      p <- c(p, -1)
    }else{
      while(status == "in progress"){
        neighs <- neighbors(net, v = src)
        closest.to.target <- neighs[which.min(hyperbolic_dist(polar[target[i], ], polar[neighs, ]))]
        while(closest.to.target != target[i] & closest.to.target %in% faulty){
          neighs <- neighs[neighs != closest.to.target]
          closest.to.target <- neighs[which.min(hyperbolic_dist(polar[target[i], ], polar[neighs, ]))]
        }
        if(closest.to.target == target[i]){
          p <- c(p, closest.to.target)
          status <- "done"
        }else if(closest.to.target == previous){
          p <- c(p, -1)
          status <- "dropped"
        }else{
          previous <- src
          src <- closest.to.target
          p <- c(p, closest.to.target)
        }
      }
    }
    paths[[i]] <- p
  }
  return(paths)
}


#' Check whether short hyperbolic distances are indicative of links
#'
#' Given a complex network and the hyperbolic coordinates of its nodes, compute the fraction of existing edges between all nodes whose pairwise
#' distances are within the considered bin.
#' 
#' @param net igraph; A complex network.
#' @param polar data frame; Polar coordinates (r, theta) of all nodes in the network.
#' @param bins numeric; The number of distance bins to consider (default = 10).
#' 
#' @return A data frame with the two following elements:
#' \item{dist}{The considered distance bins.}
#' \item{prob}{the connection probabilities within each bin.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Papadopoulos, F. et al. (2012) Popularity versus similarity in growing networks. \emph{Nature} 489(7417):537-40.
#' 
#' @examples 
#' # Generate an artificial network with the PS model
#' net <- ps_model(500, 6, 2.5, 0)
#' 
#' # Get the connection probability curve and plot it
#' conn <- get_conn_probs(net$network, net$polar, 15)
#' 
#' plot(conn$dist, conn$prob, pch = 16, xlab = "Hyperbolic distance", ylab = "Connection probability")
#' 
#' @export
#' @import igraph
#' 
get_conn_probs <- function(net, polar, bins = 10){
  
  # The network's adjacency matrix serves as the reference for edge existence
  ref <- as_adjacency_matrix(net, type = "upper", sparse = F)
  
  # Compute the matrix of pairwise hyperbolic distances between nodes
  # WARNING: This can be quite computationally intensive!!!
  dist <- sapply(seq(nrow(polar)), function(i) hyperbolic_dist(polar[i, ], polar))
  
  # Get rid of useless data
  rm(net, polar)
  
  max.dist <- max(dist)
  dist[lower.tri(dist, diag = T)] <- Inf
  ref[lower.tri(ref, diag = T)] <- 0
  
  steps <- seq(0, max.dist, length.out = bins)

  res <- data.frame(dist = steps, prob = vector("numeric", length = bins), stringsAsFactors = F)
  
  for(i in 1:(bins - 1)){
    res$prob[i] <- sum(ref[(dist >= steps[i]) & (dist < steps[i + 1])])/sum((dist >= steps[i]) & (dist < steps[i + 1]))
  }
  
  # Compute the last probability
  dist[lower.tri(dist, diag = T)] <- 0
  res$prob[bins] <- sum(ref[dist >= steps[bins]])/sum(dist >= steps[bins])
  
  return(res)
}

#' Connection probabilities in the PS model
#'
#' Given a a set of network properties, compute connection probabilities according to 
#' the Popularity-Similarity model [Papadopoulos et al. 2012, Nature 489(7417):537-40].
#' 
#' @param bins vector; The hyperbolic distance bins or steps in which theoretical connection probabilities are to be computed.
#' @param N integer; Number of network nodes.
#' @param avg.k numeric; Network's average node degree.
#' @param gma numeric; The network's scaling exponent gamma.
#' @param Temp numeric; The network temperature.
#' 
#' @return A data frame with the two following elements:
#' \item{dist}{The considered distance bins.}
#' \item{prob}{the connection probabilities within each bin.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Papadopoulos, F. et al. (2012) Popularity versus similarity in growing networks. \emph{Nature} 489(7417):537-40.
#' 
#' @examples 
#' # Generate an artificial network with the PS model
#' net <- ps_model(500, 6, 2.5, 0.1)
#' 
#' # Get the real and theoretical connection probability curves and plot them
#' conn <- get_conn_probs(net$network, net$polar, 15)
#' theo <- get_theoretical_conn_probs(conn$dist, vcount(net$network), 6, 2.5, 0.1)
#' 
#' plot(conn$dist, conn$prob, pch = 16, 
#'      xlab = "Hyperbolic distance", ylab = "Connection probability")
#' points(theo$dist, theo$prob, pch = 16, 
#'      xlab = "Hyperbolic distance", ylab = "Connection probability", col = "red")
#' legend("topright", c("Real", "Theory"), pch = c(16, 16), col = c("black", "red"))
#' 
#' @export
#' @import igraph
#' 
get_theoretical_conn_probs <- function(bins, N, avg.k, gma, Temp){
  
  # Compute beta and m
  beta <- 1 / (gma - 1) # Parameter controlling popularity fading
  m <- round(avg.k/2) # Parameter controlling average node degree
  
  # Determine the radius of the hyperbolic disc containing the network
  if(beta == 1 & Temp == 0){
    R <- 2*log(N) - 2*log(2*log(N)*(2*Temp)/(2*m*sin(Temp*pi)))
  }else if(beta == 1 & Temp > 0){
    R <- 2*log(N) - 2*log(2*log(N)*(2)/(2*m*pi))
  }else if(beta < 1 & Temp == 0){
    R <- 2*log(N) - 2*log((2*(1 - exp(-0.5*(1 - beta)*2*log(N))))/(pi*m*(1 - beta))) #SI, Eq. 10
  }else{
    R <- 2*log(N) - 2*log((2*Temp*(1 - exp(-0.5*(1 - beta)*2*log(N))))/(sin(Temp*pi)*m*(1 - beta))) #SI, Eq. 26
  }
  
  # Connection probability determination
  if(Temp == 0){
    prob <- (sign(R - bins) + 1)/2
  }else{
    prob <- 1 / (1 + exp((bins - R)/(2*Temp)))
  }
  res <- data.frame(dist = bins, prob = prob)
  
  return(res)
}