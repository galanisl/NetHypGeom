#'  Flights between the 500 busiest commercial airports in the US
#'
#' The airport network corresponds to the connections between the 500 busiest
#' commercial airports in the United States. Two airports are linked if there was a flight
#' scheduled between them in 2002.
#' The network's scaling exponent and inferred temperature are 2.01 and 0.15, respectively.
#'
#' @docType data
#'
#' @format An \code{igraph} object representing the airport network and a numeric value with its inferred temperature.
#'
#' @keywords datasets
#'
#' @references Colizza V, Pastor-Satorras R, Vespignani A (2007) Reaction-diffusion processes 
#' and metapopulation models in heterogeneous networks. Nature Physics 3:276-282. doi:10.1038/nphys560
#' @references Alanis-Lobato, G., Mier, P. & Andrade-Navarro, M. (2016) Manifold learning and 
#' maximum likelihood estimation for hyperbolic network embedding. Applied Network Science 1(10)
#'
#' @source \href{http://opsahl.co.uk/tnet/datasets/USairport500.txt}{Tore Opsahl's website}
#'
#' @examples
#' # Map the airport network to hyperbolic space using LaBNE+HM
#' coords <- labne_hm(net = air, gma = 2.01, Temp = 0.15, w = pi/12)
#' 
#' # Visually explore the resulting hyperbolic mapping
#' plot_hyperbolic_net(network = air, nodes = coords$polar, 
#'                     node.colour = coords$polar$theta)
#' 
"air"