#include <vector>
#include <cmath>
#include <string>

class Node {
  
private:
  std::string name;
  double radius, radius_init, R_i;
  double theta;
  double r_sinh, r_cosh;
  double rinit_sinh, rinit_cosh;
  int is_here;
  
  std::vector<Node*> * adjNodeList;
  
public:
  Node(){};
  Node(std::string id){ name = id; is_here=0; adjNodeList = new std::vector<Node*>(); }
  ~Node(){ delete adjNodeList; }
  
  int getK() const { return adjNodeList->size(); };
  
  void setAngle(double angle){ theta = angle; }
  double getAngle() const { return theta; }
  
  void setRadius(double r, double zeta){
    radius = r; r_sinh = sinh(zeta * r); r_cosh = cosh(zeta * r);
  }
  double getRadius() { return radius; }
  double getRsinh() { return r_sinh; }
  double getRcosh() { return r_cosh; }
  
  void setInitRadius(double r, double zeta){
    radius_init = r; rinit_sinh = sinh(zeta * r); rinit_cosh = cosh(zeta * r);
  }
  double getInitRadius() const { return radius_init; }
  double getInitsinh() const { return rinit_sinh; }
  double getInitcosh() const { return rinit_cosh; }
  
  void setR(double R){ R_i = R; };
  double getR() { return R_i; }
  
  void set_is_here() {is_here = 1;}
  int get_is_here() {return is_here;}
  
  std::string getName() const { return name; }
  
  void addAdjNode(Node *adj) { adjNodeList->push_back(adj); }
  
  std::vector<Node*> * getAdjNodeList() const { return adjNodeList; }
  
  bool isNeighbor(Node *x) const {
    if(std::find(adjNodeList->begin(), adjNodeList->end(), x) != adjNodeList->end())
      return true;
    else
      return false;
  }
};
