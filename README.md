NetHypGeom
================
Gregorio Alanis-Lobato, Pablo Mier and Miguel A. Andrade-Navarro

Tools for the generation, analysis and visualisation of complex networks in hyperbolic space.

Introduction
============

The network representation of many complex systems, like the Internet or the protein interactome, shows characteristics commonly present in geometric objects; scale invariance and self-similarity amongst them<sup>[1](#ref-Barabasi1999)–[4](#ref-Serrano2008)</sup>. It is then no surprise that several models, aimed at mimicking the evolution and formation of these networks, assume the existence of a hidden geometry underlying their structure and shaping their topology<sup>[5](#ref-Barthelemy2011)</sup>.

Of special interest is the so-called Popularity-Similarity (PS) model, which sustains that strong clustering and scale invariance are the result of an optimisation process involving two measures of attractiveness: node popularity and similarity between nodes<sup>[6](#ref-Papadopoulos2012)</sup>. On the one hand, popularity reflects the ability of a node to attract connections from others over time, and it is thus associated with its seniority status in the system. On the other, nodes that are similar simply tend to connect, regardless of their rank.

The PS model has a geometric interpretation in hyperbolic space, where the trade-offs that new nodes have to optimise when joining a system are abstracted by the hyperbolic distance between them and existing ones<sup>[6](#ref-Papadopoulos2012),[7](#ref-Krioukov2010)</sup>. If the PS model can generate networks that are similar to those we observe in nature and engineering, does it mean that packets travelling the Internet, signals going from receptors to transcription factors in the cell or messages between people in social networks traverse the hyperbolic geometry underlying each of these systems? To answer this question, we need a means to map them to hyperbolic space.

In 2015, Papadopoulos and colleagues introduced HyperMap, a Maximum Likelihood Estimation (MLE) approach, in which the space of PS models with the same structural properties as the network of interest is explored, in search for the one that better fits its topology<sup>[8](#ref-Papadopoulos2015)</sup>. This search is very accurate, albeit computationally demanding.

Inspired by the well-established field of non-linear dimensionality reduction in Machine Learning<sup>[9](#ref-Cayton2005)</sup>, we recently put forward the Laplacian-based Network Embedding or LaBNE<sup>[10](#ref-Alanis2016a)</sup>. LaBNE extremely fast, but highly depends on topological information to carry out good embeddings. In addition, LaBNE’s aim to map connected nodes as close as possible to each other in the embedding space, disregarding that disconnected nodes should be far apart<sup>[11](#ref-Shaw2009)</sup>, can lead to inaccuracies when associating short hyperbolic distances with connections between nodes.

Given the drawbacks and limitations of both HyperMap and LaBNE, we aimed at combining them to improve LaBNE's accuracy and reduce HyperMap's execution times. LaBNE+HyperMap (LaBNE+HM), our proposed approach, uses LaBNE to quickly draft a first geometric configuration of a network of interest. This draft is then passed on to HyperMap, which refines the embedding and produces the final mapping to hyperbolic space<sup>[12](#ref-Alanis2016b)</sup>. LaBNE+HM profits from LaBNE's fast embeddings and significantly reduces the search space of HyperMap.

With `NetHypGeom` it is possible to generate artificial networks using the above-described PS model, visualise these networks in hyperbolic space and explore their topological and geometric characteristics. In addition, it is possible to embed artifical and real networks to hyperbolic space using LaBNE+HM, LaBNE or HyperMap and study how good these mappings are. LaBNE and the manifold learning part of LaBNE+HM were fully implemented by us in R. HyperMap and the MLE part of LaBNE+HM were implemented by the [DK Lab](https://bitbucket.org/dk-lab/2015_code_hypermap) in C++. We wrote an R wrapper that consumes this code when needed.

Installation
============

1.  Install the `devtools` package from CRAN if you haven't done so:

``` r
install.packages("devtools")
```

1.  Load the `devtools` package:

``` r
library("devtools")
```

1.  Install `NetHypGeom` using the `install_github` function:

``` r
install_github("galanisl/NetHypGeom")
```

Usage
=====

To start using `NetHypGeom`, load the package:

``` r
library("NetHypGeom")
```

Let's generate an artificial network with the PS model, visualise its node degree distribution and compare its clustering coefficient with that of random networks. The network will have *N* = 500 nodes, an average node degree $\\bar{k} = 10$, a target scaling exponent *γ* = 2.3 and a strong clustering coefficient, i.e. a low temperature *T* = 0.15<sup>[6](#ref-Papadopoulos2012)</sup>:

``` r
net <- ps_model(N = 500, avg.k = 10, gma = 2.3, Temp = 0.15)
plot_degree_distr(net$network)
compare_clustering_to_er(net$network, 100, "PS network")
```

Let's now access the coordinates of nodes to visualise the network in hyperbolic space, using the node angular (similarity) dimension as colour:

``` r
plot_hyperbolic_net(network = net$network, nodes = net$polar, node.colour = net$polar$theta)
```

To investigate whether short hyperbolic distances do correspond to the presence of links, we can explore connection probabilities over the range of distances between network nodes, using a specified number of distance bins:

``` r
conn <- get_conn_probs(net = net$network, polar = net$polar, bins = 15)
plot(conn$dist, conn$prob, pch = 16, xlab = "Hyperbolic distance", ylab = "Connection probability")
```

One of the big advantages of revealing the geometry underlying a complex network is that it enables the analysis of its navigation efficiency. An important function of complex systems is the routing of information or signals (that we refer to as packets here) without global knowledge of the network topology, avoiding loss of the packet and following short paths<sup>[13](#ref-Boguna2009),[14](#ref-Papadopoulos2010)</sup>. With `NetHypGeom` we can check if it is possible to efficiently send packets from a source node to a target one using greedy routing. This means that the source node ships a packet to the direct neighbour that is hyperbolically closest to the target, the recipient neighbour does the same with its direct neighbours and so on, until the packet reaches the target. If, in the delivery process, a neighbour sends the packet to the previously visited node, i.e. it falls into a loop, the packet is dropped and the delivery is flagged as unsuccessful.

``` r
# Form vectors of random non-redundant source-target pairs
N <- vcount(net$network)
st <- 1000

# We subtract 1, because the formulae to go from linear upper 
# diagonal indexing to (i,j) are zero-based
k <- sample(N*(N-1)/2, st) - 1
sources <- (N - 2 - floor(sqrt(-8*k + 4*N*(N-1)-7)/2.0 - 0.5))
targets <- (k + sources + 1 - N*(N-1)/2 + (N-sources)*((N-sources)-1)/2)

# Back to 1-based indexing
sources <- sources + 1
targets <- targets + 1

# Analyse the network's navigability
hop.stretch <- greedy_route_packets(net$network, net$polar, sources, targets)

# Compute the fraction of succesfully delivered packets
sum(hop.stretch)/st
```

Finally, let's take advantage of the fact that we know the hyperbolic coordinates of nodes in our artificial network to compare the mapping to hyperbolic space using LaBNE+HM and HyperMap:

``` r
# To embed the network using HyperMap, we set LaBNE+HM's window to 2*pi
hm <- labne_hm(net = net$network, gma = 2.3, Temp = 0.15, k.speedup = 10, w = 2*pi)

# To embed with LaBNE+HM, we reduce HyperMap's search space from 2*pi 
# to a small window of 15 degrees around LaBNE's angles
lh <- labne_hm(net = net$network, gma = 2.3, Temp = 0.15, k.speedup = 10, w = pi/12)

# Comparison between real and HyperMap-inferred angles and real and LaBNE+HM-inferred angles
plot(net$polar$theta, hm$polar$theta, pch = 16, 
     xlab = "Real angles", ylab = "Inferred angles", main = "HyperMap")
plot(net$polar$theta, lh$polar$theta, pch = 16, 
     xlab = "Real angles", ylab = "Inferred angles", main = "LaBNE+HM")
```

How to cite
===========

If you find this package useful, please cite the following publications:

-   Alanis-Lobato, G., Mier, P. & Andrade-Navarro, M. Manifold learning and maximum likelihood estimation for hyperbolic network embedding. *Applied Network Science* **1**(10) (2016) [See paper](http://appliednetsci.springeropen.com/articles/10.1007/s41109-016-0013-0)
-   Alanis-Lobato, G., Mier, P. & Andrade-Navarro, M. Efficient embedding of complex networks to hyperbolic space via their Laplacian. *Scientific Reports* **6**, 30108 (2016) [See paper](http://www.nature.com/articles/srep30108)
-   Papadopoulos, F., Aldecoa, R. & Krioukov, D. Network geometry inference using common neighbors. *Physical Review E* **92**, 022807 (2015) [See paper](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.022807)

References
==========

1. Barabási, A.-L. & Albert, R. Emergence of scaling in random networks. *Science* **286,** 509–512 (1999).

2. Song, C., Havlin, S. & Makse, H. A. Origins of fractality in the growth of complex networks. *Nat. Phys.* **2,** 275–281 (2006).

3. Goh, K.-I., Salvi, G., Kahng, B. & Kim, D. Skeleton and fractal scaling in complex networks. *Phys. Rev. Lett.* **96,** 018701 (2006).

4. Serrano, M. Á., Krioukov, D. & Boguñá, M. Self-similarity of complex networks and hidden metric spaces. *Phys. Rev. Lett.* **100,** 078701 (2008).

5. Barthélemy, M. Spatial networks. *Phys. Rep.* **499,** 1–101 (2011).

6. Papadopoulos, F., Kitsak, M., Serrano, M. Á., Boguñá, M. & Krioukov, D. Popularity versus similarity in growing networks. *Nature* **489,** 537–540 (2012).

7. Krioukov, D., Papadopoulos, F., Kitsak, M., Vahdat, A. & Boguñá, M. Hyperbolic geometry of complex networks. *Phys. Rev. E* **82,** 036106 (2010).

8. Papadopoulos, F., Aldecoa, R. & Krioukov, D. Network geometry inference using common neighbors. *Phys. Rev. E* **92,** 022807 (2015).

9. Cayton, L. Algorithms for manifold learning. *UCSD tech report* **CS2008-0923,** 1–17 (2005).

10. Alanis-Lobato, G., Mier, P. & Andrade-Navarro, M. A. Efficient embedding of complex networks to hyperbolic space via their Laplacian. *Sci. Rep.* **6,** 30108 (2016).

11. Shaw, B. & Jebara, T. Structure preserving embedding. in *Proceedings of the 26th Annual International Conference on Machine Learning* 937–944 (ACM, 2009).

12. Alanis-Lobato, G., Mier, P. & Andrade-Navarro, M. A. Manifold learning and maximum likelihood estimation for hyperbolic network embedding. *Applied Network Science* **1,** 10 (2016).

13. Boguñá, M., Krioukov, D. & Claffy, K. C. Navigability of complex networks. *Nat. Phys.* **5,** 74–80 (2009).

14. Papadopoulos, F., Krioukov, D., Boguñá, M. & Vahdat, A. Greedy forwarding in dynamic scale-free networks embedded in hyperbolic metric spaces. in *INFOCOM, 2010 proceedings IEEE* 1–9 (2010).
