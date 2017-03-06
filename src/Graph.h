#include <cstdio>
#include "Node.h"
#include <Rcpp.h>

using namespace Rcpp;

class Graph {
  
private:
  double kbar;
  double m;
  double L;
  std::vector<Node*> * nodeList;
  
  // function to sort by degree
  static bool compByK(Node *a, Node *b);
  
public:
  Graph();
  ~Graph();
  
  // Set parameters
  void setM(double value){ m = value; } // m specified by the user 
  void setL(double value){ L = value; } // L specified by the user 
  
  // Get parameters
  double getN() const { return nodeList->size(); }
  double getM() const { return m; }
  double getL() const { return L; }
  double getKbar() const { return kbar; }
  std::vector<Node*> * getNodeList() const { return nodeList; }
  
  // Read and write graph functions
  bool readGraph(DataFrame net);
  bool readLaBNEAngles(DataFrame theta);
  DataFrame getCoords();
  void addNewNode(Node *node) { nodeList->push_back(node); }
  Node * findNode(const std::string& name);
};
