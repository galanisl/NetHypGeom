#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>

#include "Graph.h"

using namespace std;

Graph::Graph(){
  m = -1;
  L = -1;
  nodeList = new std::vector<Node*>;
}

Graph::~Graph(){
  std::vector<Node*>::iterator it = nodeList->begin();
  for(; it != nodeList->end(); ++it)
    delete *it;
  
  delete nodeList;
}


bool Graph::compByK(Node *a, Node *b){
  return (int)a->getK() > (int)b->getK();
}


Node * Graph::findNode(const std::string& name) {
  for (int i = 0; i < nodeList->size(); i++) {
    std::string s = nodeList->at(i)->getName(); 
    if (s.compare(name) == 0)
      return nodeList->at(i);
  }
  return NULL;
}

bool Graph::readLaBNEAngles(DataFrame theta){
  CharacterVector nodeID = theta["id"];
  NumericVector tht = theta["theta"];
  std::string id;
  double ang;
  for(int i = 0; i < theta.nrows(); i++){
    id = nodeID[i];
    ang = tht[i];
    Node *u = findNode(id);
    if(u != NULL){
      u->setAngle(ang);
    }
  }
  return true;
}

bool Graph::readGraph(DataFrame net){
  
  CharacterVector vertex1 = net[0];
  CharacterVector vertex2 = net[1];
  std::string v1, v2;
  int i = 0;
  for(i = 0; i < net.nrows(); i++){
    v1 = vertex1[i];
    Node *u = findNode(v1);
    if (u == NULL) {
      u = new Node(v1);
      addNewNode(u);
    }
    v2 = vertex2[i];
    Node *v = findNode(v2);
    if (v == NULL) {
      v = new Node(v2);
      addNewNode(v);
    }
    
    u->addAdjNode(v);
    v->addAdjNode(u);
  }
  
  // sort by degree
  std::sort(nodeList->begin(), nodeList->end(), compByK);
  
  // Calculate L and m
  int totalK = 0;
  for(i = 0; i < nodeList->size(); ++i)
    totalK += nodeList->at(i)->getK();
  kbar = (double)totalK / nodeList->size();
  
  if(m == -1) // if not specified by the user, set it equal to the minimum node degree.
    m = nodeList->at(nodeList->size()-1)->getK();
  if(m < 0) {
    m = 0;
  }
  
  if(L == -1) // if not specified by the user.
    L = (kbar - 2.0 * m) / 2.0;
  if(L < 0) {
    L = 0;
  }
  
  return true;
}

DataFrame Graph::getCoords(){
  std::vector<std::string> id;
  std::vector<double> r;
  std::vector<double> theta;
  for (int i = 0; i < nodeList->size(); ++i) {
    Node *u = nodeList->at(i);
    id.push_back(u->getName());
    r.push_back(u->getRadius());
    theta.push_back(u->getAngle());
  }
  DataFrame coords = DataFrame::create(Named("id")=id, Named("r")=r, Named("theta")=theta);
  return(coords);
}