# NetHypGeom
Tools for the generation, analysis and visualisation of complex networks in hyperbolic space.

## Introduction

The network representation of many complex systems, like the Internet or the protein interactome, shows characteristics commonly present in geometric objects; scale invariance and self-similarity amongst them (Barabási and Albert 1999; Song et al. 2006; Goh et al. 2006; Serrano et al. 2008). It is then no surprise that several models, aimed at mimicking the evolution and formation of these networks, assume the existence of a hidden geometry underlying their structure and shaping their topology (Barthélemy 2011).

Of special interest is the so-called Popularity-Similarity (PS) model, which sustains that strong clustering and scale-free node degree distributions are the result of an optimisation process involving two measures of attractiveness: node popularity and similarity between nodes (Papadopoulos et al. 2012). On the one hand, popularity reflects the ability of a node to attract connections from other nodes over time, and it is thus associated with its seniority status in the system. On the other, nodes that are similar simply tend to connect, regardless of their rank.

The PS model has a geometric interpretation in hyperbolic space, where the trade-offs that new nodes have to optimise when joining a system are abstracted by the hyperbolic distance between them and existing ones (Krioukov et al. 2010; Papadopoulos et al. 2012). If the PS model can generate networks that are similar to those we observe in nature and engineering (Krioukov et al. 2010; Papadopoulos et al. 2012), does it mean that packets travelling the Internet, signals going from receptors to transcription factors in the cell or messages between people in social networks traverse the hyperbolic geometry underlying each of these systems? To answer this question, we need a means to map them to hyperbolic space.

In 2015, Papadopoulos and colleagues introduced HyperMap, a Maximum Likelihood Estimation (MLE) approach, in which the space of PS models with the same structural properties as the network of interest is explored, in search for the one that better fits its topology (see Methods and (Papadopoulos F et al. 2015b; Papadopoulos et al. 2015a) for more details). This search is very accurate, albeit computationally demanding.

Inspired by the well-established field of non-linear dimensionality reduction in Machine Learning (Cayton 2005), we recently put forward the Laplacian-based Network Embedding or LaBNE (Alanis-Lobato et al. 2016). LaBNE extremely fast (see Fig. 1a), but highly depends on topological information to carry out good embeddings.  In addition, LaBNE’s aim to map connected nodes as close as possible to each other in the embedding space, disregarding that disconnected nodes should be far apart (Shaw and Jebara 2009), can lead to inaccuracies when associating short hyperbolic distances with connections between nodes.

Given the drawbacks and limitations of both HyperMap and LaBNE, we aimed at combining them to improve LaBNE’s accuracy and reduce HyperMap’s execution times. LaBNE+HyperMap (LaBNE+HM), our proposed approach, uses LaBNE to quickly draft a first geometric configuration of a network of interest in H2 . This draft is then passed on to HyperMap, which refines the embedding and produces the final mapping to hyperbolic space. LaBNE+HM profits from LaBNE’s fast embeddings and significantly reduces the search space of HyperMap.

## Installation

1. Install the `devtools` package from CRAN if you haven't done so:

```r
install.packages("devtools")
```

2. Load the `devtools` package:

```r
library("devtools")
```

3. Install `NetHypGeom` using the `install_github` function:

```r
install_github("galanisl/NetHypGeom")
```

## Usage

## How to cite

If you find this package useful, please cite the following publications:

- Alanis-Lobato, G., Mier, P. & Andrade-Navarro, M. (2016) Manifold learning and maximum likelihood estimation for hyperbolic network embedding. *Applied Network Science* 1(10) [See paper](http://appliednetsci.springeropen.com/articles/10.1007/s41109-016-0013-0)
- Alanis-Lobato, G., Mier, P. & Andrade-Navarro, M. (2016) Efficient embedding of complex networks to hyperbolic space via their Laplacian. *Scientific Reports* 6, 30108 [See paper](http://www.nature.com/articles/srep30108)
- Papadopoulos, F., Aldecoa, R. & Krioukov, D. (2015) Network geometry inference using common neighbors. *Physical Review E* 92, 022807 [See paper](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.022807)

