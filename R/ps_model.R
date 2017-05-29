
#' Generate a Popularity-Similarity network
#' 
#' Generates an artificial network following the Popularity-Similarity model proposed in [Papadopoulos et al. 2012, Nature 489(7417):537-40].
#' This implementation only considers the hyperbolic space curvature \eqn{K = -(1^2)}.
#' 
#' @param N integer; The number of desired nodes for the network (default = 500).
#' @param avg.k numeric; The target average node degree (default = 10).
#' @param gma numeric; The target scaling exponent of the network's node degree distribution (default = 2).
#' @param Temp numeric; The network temperature, controlling the target clustering coefficient (default = 0).
#' 
#' @return List with the two following elements:
#' \item{network}{igraph object representation of the input network.}
#' \item{polar}{Data frame containing elements \code{r} and \code{theta}, the radial and angular coordinates of the network nodes in hyperbolic space.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Papadopoulos, F. et al. (2012) Popularity versus similarity in growing networks. \emph{Nature} 489(7417):537-40.
#' 
#' @examples
#' net.default <- ps_model()
#' net <- ps_model(100, 6, 2.5, 0.2)
#' 
#' @export
#' @import igraph
#' 
ps_model <- function(N = 500, avg.k = 10, gma = 2, Temp = 0){
  
  # Make sure that all parameters are valid
  check_ps_parameters(N, avg.k, gma, Temp)
  
  beta <- 1 / (gma - 1) # Determine the value of beta, a parameter controlling popularity fading
  m <- round(avg.k/2)
  
  # Initialise the node data frame and add the first node to the network
  nodes <- data.frame(r = vector("numeric", length = N), theta = vector("numeric", length = N))
  nodes$r[1] <- 0
  nodes$theta[1] <- stats::runif(1, min = 0, max = 2*pi)
  
  # Initialise the network adjacency matrix
  net <- matrix(0, nrow = N, ncol = N)
  diag(net) <- 1 #Trick to avoid a node choosing itself as interactor in the process
  
  # Dummy variable
  m_count <- 0
  
  # Add the rest of the nodes and connections to the network
  for(t in 2:N){
    
    # Move all nodes to their new radial coordinates to simulate popularity fading
    nodes$r[1:(t - 1)] <- beta * 2 * log(1:(t - 1)) + (1 - beta) * 2 * log(t)
    
    # New node is added to the network and acquires polar coordinates
    nodes$r[t] <- 2*log(t)
    nodes$theta[t] <- stats::runif(1, min = 0, max = 2*pi)
    
    # If beta = 1, popularity fading is not simulated
    if(beta == 1 & Temp == 0){
      R.t <- 2*log(t) - 2*log(2*log(t)*(2*Temp)/(2*m*sin(Temp*pi)))
    }else if(beta == 1 & Temp > 0){
      R.t <- 2*log(t) - 2*log(2*log(t)*(2)/(2*m*pi))
    }else if(beta < 1 & Temp == 0){
      R.t <- 2*log(t) - 2*log((2*(1 - exp(-0.5*(1 - beta)*2*log(t))))/(pi*m*(1 - beta))) #SI, Eq. 10
    }else{
      R.t <- 2*log(t) - 2*log((2*Temp*(1 - exp(-0.5*(1 - beta)*2*log(t))))/(sin(Temp*pi)*m*(1 - beta))) #SI, Eq. 26
    }
    
    # Hyperbolic distance between new node and all existing nodes
    d <- hyperbolic_dist(nodes[t, ], nodes[1:(t - 1), ])
    
    if(Temp > 0 & t > 1 & (t - 1) > m){
      # Compute the probability that the new node t connects to all existing nodes
      p = 1 / (1 + exp((d - R.t)/(2*Temp)))
      while(m_count < m){
        if(length(which(net[t, 1:(t - 1)] == 0)) <= m){
          net[t, which(net[t, 1:(t - 1)] == 0)] <- 1
          m_count <- m
        }else if(sum(p[which(net[t, 1:(t - 1)] == 0)] < 0.1) == length(which(net[t, 1:(t - 1)] == 0))){
          s <- sort(d[which(net[t, 1:(t - 1)] == 0)], index.return = T)
          if(length(d[which(net[t, 1:(t - 1)] == 0)]) < m){
            net[t, s$ix] <- 1
          }else{
            net[t, s$ix[1:m]] <- 1
          }
          m_count <- m
        }
        else{
          # Sample m nodes randomly that are not partners of the new node yet
          rnd_nodes <- sample(which(net[t, 1:(t - 1)] == 0), m - m_count)
          # Connect to these m nodes with probability p
          connect.to <- stats::runif(m - m_count) <= p[rnd_nodes]
          net[t, rnd_nodes[connect.to]] <- 1
          m_count <- m_count + sum(connect.to)
        }
      }
      m_count <- 0
    }else{
      # Since Temp = 0, simply connect to the m hyperbolically closest nodes
      s <- sort(d, index.return = T)
      if(length(d) < m){
        net[t, s$ix] <- 1
      }else{
        net[t, s$ix[1:m]] <- 1
      }
    }
  }
  diag(net) <- 0
  net <- net + t(net)
  net <- graph_from_adjacency_matrix(net, mode = "undirected", diag = FALSE)
  return(list(network = net, polar = nodes))
}

#' Check PS model parameters
#' 
#' Function that makes sure that all parameters passed to ps.model are valid.
#' 
#' @param N The number network nodes.
#' @param avg.k The target average node degree.
#' @param gma The network's scaling exponent.
#' @param Temp The network temperature.
#' 
#' @return None
#' 
check_ps_parameters <- function(N, avg.k, gma, Temp){
  
  ### N (Number of network nodes)
  if(is.na(N) | N < 2){
    stop("The number of network nodes N was not specified or is not a valid integer >= 2.")
  }
  
  ### AVG.K (Average node degree)
  if(is.na(avg.k) | avg.k < 2){
    stop("The target average node degree avg.k was not specified or is not a valid number >= 2.")
  }
  
  ### GAMMA (Network's scaling exponent)
  if(is.na(gma) | gma < 2 | gma > 3){
    stop("The target network scaling exponent gma was not specified or is outside the valid range [2, 3].")
  }
  
  ### TEMPERATURE (Network temperature)
  if(is.na(Temp)){
    stop("You have to specify a temperature value in the range [0, Inf).")
  }else if(Temp < 0){
    stop("The specified temperature value is outside the valid range [0, Inf). Please specify a valid value.")
  }
}
