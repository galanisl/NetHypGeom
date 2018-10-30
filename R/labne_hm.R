
#' Laplacian-Based Network Embedding+HyperMap (LaBNE+HM)
#' 
#' LaBNE+HM is a method for the fast yet accurate embedding of complex networks to the native representation of the two-dimensional hyperbolic plane.
#' 
#' @param net igraph, data frame or path to tab-separated file; The complex network to be embedded into hyperbolic space.
#' @param gma numeric; The network's scaling exponent. If not specified, it is automatically computed.
#' @param Temp numeric; The network's temperature (low temperature for strongly clustered networks and vice versa). If not specified, it is set to 0.1.
#' @param k.speedup integer; A speedup heuristic will be applied to nodes of degree < k.speedup. Its default value is 10. If set to 0, no heuristic is applied.
#' @param m.in numeric; The expected initial node degree, i.e. the average number of link stubs with which a node joins the network.  Parameter m.in can be obtained from historical data of the evolution of the network. If this data is not available and m.in is not specified, it is set to the minimum observed node degree in the network.
#' @param L.in numeric; The internal link formation rate, i.e. the average number of links formed between existing network nodes. If not specified, it is set to L = (kbar-2*m)/2, where kbar is the average node degree of the network.
#' @param w numeric; The window considered by HyperMap to refine angles found by LaBNE. If not specified, it is set to "auto", which means that w is set to 2*pi*Temp^2.
#' 
#' @return List with the three following elements:
#' \item{network}{igraph object representation of the input network.}
#' \item{polar}{Data frame containing elements \code{r} and \code{theta}, the inferred radial and angular coordinates of the network nodes in hyperbolic space.}
#' \item{cartesian}{Data frame containing elements \code{x} and \code{y}, the inferred cartesian coordinates of the network nodes in hyperbolic space. It is useful to have these coordinates for the direct visualisation of the embedding using standard plotting function in R.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Alanis-Lobato, G., Mier, P. and Andrade-Navarro, M. (2016) Manifold learning and maximum likelihood estimation for hyperbolic network embedding. \emph{Applied Network Science} 1(10).
#' @references Alanis-Lobato, G., Mier, P. and Andrade-Navarro, M. (2016) Efficient embedding of complex networks to hyperbolic space via their Laplacian. \emph{Scientific Reports} 6, 30108.
#' 
#' @examples
#' # Generate an artificial network with the PS model, 
#' # such that the hyperbolic coordinates of its nodes are known
#' net <- ps_model(N = 500, avg.k = 10, gma = 2.5, Temp = 0.15)
#' 
#' # Map the network to hyperbolic space using LaBNE+HM and specifying all parameters.
#' coords <- labne_hm(net = net$network, gma = 2.5, Temp = 0.15, 
#'                    k.speedup = 10, m.in = 5, L.in = 0, w = pi/12)
#' 
#' # Visually explore the resulting hyperbolic mapping
#' plot_hyperbolic_net(network = net$network, 
#'                     nodes = coords$polar, node.colour = net$polar$theta)
#' 
#' @useDynLib NetHypGeom, .registration = TRUE
#' 
#' @export
#' @import igraph
#' @import Rcpp
#' @importFrom RSpectra eigs_sym
#' 
labne_hm <- function(net, gma = NA, Temp = 0.1, k.speedup = 10, m.in = NA, L.in = NA, w = "auto"){

  options(stringsAsFactors = F)
    
  # Make sure that all parameters are valid
  params <- check_lhm_parameters(net, gma, Temp, k.speedup, m.in, L.in, w)
  
  N <- vcount(params$net) #The number of nodes
  bta <- 1 / (params$gma - 1) #Determine the value of beta, a parameter controlling popularity fading
  
  # Perform Laplacian Eigenmaps (projection over the two eigenvectors associated with the smallest non-zero eigenvalues)
  # The algorithm of RSpectra is good at finding large eigenvalues rather than small ones, so setting which = "SM" (smallest magnitude) 
  # tends to require much more iterations or even may not converge. The recommended way to retrieve the smallest eigenvalues is to
  # calculate the largest eigenvalues of A^-1, whose reciprocals are exactly the smallest eigenvalues of A. 
  # RSpectra implements this through parameter sigma.
  # See https://cran.r-project.org/web/packages/RSpectra/vignettes/introduction.html for more details.
  leigs <- eigs_sym(laplacian_matrix(params$net, weights = NA), k = 3, which = "LM", sigma = 0)
  
  # Take the eigenvectors corresponding to the smallest eigenvalues to determine the nodes' Cartesian coordinates
  cart <- data.frame(id = V(params$net)$name, x = leigs$vectors[, 2], y = leigs$vectors[, 1])
  
  # Now determine LaBNE-based node angles
  theta = atan2(cart$y, cart$x)
  theta[theta < 0] <- theta[theta < 0] + 2*pi
  
  # If w is 0, no refinement is needed
  if(w == 0){
    radial <- 2*bta*log(1:N) + 2*(1 - bta)*log(N)
    radial[sort(degree(params$net), decreasing = T, index.return = T)$ix] <- radial
    polar <- data.frame(id = V(params$net)$name, r = radial, theta = theta)
  }else{
    # Refine LaBNE-based angles and determine node hyperbolic coordinates with HyperMap
    theta <- data.frame(id = V(params$net)$name, theta = theta)
    polar <- hypermap(as_data_frame(params$net, what = "edges"), params$gma, params$Temp, params$k.speedup, params$m.in, params$L.in, params$w, theta)
    # Ensure that all angles are in [0, 2pi]
    polar$theta[polar$theta > 2*pi] <- polar$theta[polar$theta > 2*pi] - 2*pi
    polar$theta[polar$theta < 0] <- polar$theta[polar$theta < 0] + 2*pi
    cart$id <- polar$id
  }
  # Readjust the cartesian coordinates
  cart$x <- polar$r*cos(polar$theta)
  cart$y <- polar$r*sin(polar$theta)
  
  # Ensure that the resulting data frame have the same node order as the input network
  nodes <- data.frame(id = V(params$net)$name)
  cart <- merge(nodes, cart, by = "id", sort = FALSE)
  polar <- merge(nodes, polar, by = "id", sort = FALSE)
  
  return(list(network = net, polar = polar, cartesian = cart))
}


#' Check LaBNE+HM parameters
#' 
#' Function that makes sure that all parameters passed to labne_hm are valid.
#' 
#' @param net The complex network to be embedded into hyperbolic space.
#' @param gma The network's scaling exponent.
#' @param Temp The network's temperature.
#' @param k.speedup A speedup heuristic will be applied to nodes of degree < k.speedup.
#' @param m.in The expected initial node degree.
#' @param L.in The internal link formation rate.
#' @param w The window considered by HyperMap to refine angles found by LaBNE.
#' 
#' @return A list with all the revised parameters.
#' 
check_lhm_parameters <- function(net, gma, Temp, k.speedup, m.in, L.in, w){
  
  ### THE NETWORK
  if(missing(net)){
    stop("You have to at least specify a network to embed. This can be an igraph object, a data frame or a tab-separated file.")
  }else{
    if("character" %in% class(net)){ #The network is a tab-separated file
      net <- utils::read.table(net, header = F, sep = "\t", quote = "", stringsAsFactors = F)
      if(ncol(net) < 2){
        stop("The specified file must have at least two tab-separated columns.")
      }else{
        #Consider the first two columns only and construct an igraph object
        net <- net[, 1:2]
        net <- graph_from_data_frame(net, directed = F)
        if(!is_connected(net) | !is_simple(net)){
          warning("The specified igraph object represents a disconnected network, has parallel edges or some nodes have self-interactions. 
            Only the largest connected component will be considered and/or loops/parallel edges will be removed.")
          comps <- decompose(net)
          lcc_idx <- which.max(sapply(comps, vcount))
          net <- simplify(comps[[lcc_ids]], edge.attr.comb = "min")
        }
      }
    }else if(!any(class(net) %in% c("igraph", "tbl_graph", "data.frame", "tbl_df", "tbl"))){
      stop("The specified network is neither an igraph object nor a data frame or a valid tab-separated file.")
    }else if(any(class(net) %in% c("data.frame", "tbl_df", "tbl"))){
      #Consider the first two columns only and construct an igraph object
      net <- net[, 1:2]
      net <- graph_from_data_frame(net, directed = F)
      if(!is_connected(net) | !is_simple(net)){
        warning("The specified igraph object represents a disconnected network, has parallel edges or some nodes have self-interactions. 
            Only the largest connected component will be considered and/or loops/parallel edges will be removed.")
        comps <- decompose(net)
        lcc_idx <- which.max(sapply(comps, vcount))
        net <- simplify(comps[[lcc_ids]], edge.attr.comb = "min")
      }
    }else if(any(class(net) %in% c("igraph", "tbl_graph")) & (!is_connected(net) | !is_simple(net))){
      warning("The specified igraph object represents a disconnected network, has parallel edges or some nodes have self-interactions. 
            Only the largest connected component will be considered and/or loops/parallel edges will be removed.")
      comps <- decompose(net)
      lcc_idx <- which.max(sapply(comps, vcount))
      net <- simplify(comps[[lcc_ids]], edge.attr.comb = "min")
    }
    if(is.null(V(net)$name)){
      V(net)$name <- 1:vcount(net)
    }
  }
  
  ### GAMMA (Network's scaling exponent)
  if(is.na(gma)){
    gma <- fit_power_law(degree(net))$alpha
  }
  if(gma < 2 | gma > 3){
    stop("The specified or computed value of gamma is outside the valid range [2, 3]. Please specify a valid value.")
  }
  
  ### TEMPERATURE (Network temperature)
  if(is.na(Temp)){
    stop("You have to specify a temperature value in the range [0, Inf).")
  }else if(Temp < 0){
    stop("The specified temperature value is outside the valid range [0, Inf). Please specify a valid value.")
  }else if(Temp == 0){
    Temp <- 0.0001
  }
  
  ### K_SPEEDUP (Heuristic to improve HyperMap's performance)
  if(is.na(k.speedup)){
    stop("You have to specify a k.speed value in the range [0, Inf].")
  }else if(k.speedup < 0){
    stop("The specified k.speedup value is outside the valid range [0, Inf]. Please specify a valid value.")
  }else if(k.speedup > max(degree(net))){
    warning("The specified k.speedup value is larger than the maximum node degree. The speedup heuristic will be applied to all nodes. This might result in an inaccurate embedding.")
    k.speedup <- max(degree(net))
  }
  
  ### M.IN (Expected initial node degree)
  if(is.na(m.in)){
    m.in <- min(degree(net))
  }
  if(m.in < 0 | m.in > max(degree(net))){
    stop("The specified value of m.in is outside the valid range [0, Inf]. Please specify a valid value.")
  }
  
  ### L.IN (Initial link formation rate)
  if(is.na(L.in)){
    L.in <- (mean(degree(net))-2*m.in)/2
    #If the resulting value of L.in is negative, set it to 0
    if(L.in < 0){
      L.in <- 0
    } 
  }
  if(L.in < 0){
    stop("The specified value of L.in is outside the valid range [0, Inf]. Please specify a valid value.")
  }
  
  ### W (Angle refinement window)
  if(is.na(w)){
    stop("You have to specify a w value in the range [0, 2pi].")
  }else if(class(w) == "character" & w != "auto"){
    stop("You have to specify a w value in the range [0, 2pi] or set it to 'auto'.")
  }else if(w == "auto"){
    w <- 2*pi*Temp^2
  }else if(w < 0 | w > 2*pi){
    stop("The specified value of w is outside the valid range [0, 2pi]. Please specify a valid value.")
  }
  
  parameters <- list(net = net, gma = gma, Temp = Temp, k.speedup = k.speedup, m.in = m.in, L.in = L.in, w = w)
  return(parameters)
}
