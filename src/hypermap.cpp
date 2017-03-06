#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <Rcpp.h>

#include "Graph.h"
#include "gauss_legendre.h"

using namespace std;

// Size of the intersection between two vectors
int intersection(vector<Node*>& v1, vector<Node*>& v2){
  vector<Node*> v3;
  sort(v1.begin(), v1.end());
  sort(v2.begin(), v2.end());
  set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
  return v3.size();
}

struct toIntegrate
{
  double zeta;
  double zetaOver2T;
  double theta_u;
  double theta_v;
  double ruk1_sinh;
  double ruk1_cosh;
  double ruk2_sinh;
  double ruk2_cosh;
  double rvk1_sinh;
  double rvk1_cosh;
  double rvk2_sinh;
  double rvk2_cosh;
  double R_u_k;
  double R_v_k;
  
  inline double operator()(double x) const {
    double dtheta1 = M_PI-fabs(M_PI-fabs(x-theta_v));
    double dtheta2 = M_PI-fabs(M_PI-fabs(x-theta_u));
    double x_v = (1./zeta)*acosh(rvk1_cosh*rvk2_cosh - rvk1_sinh*rvk2_sinh*cos(dtheta1));
    double x_u = (1./zeta)*acosh(ruk1_cosh*ruk2_cosh - ruk1_sinh*ruk2_sinh*cos(dtheta2));
    double first = 1. / (1. + exp(zetaOver2T * (x_v - R_v_k)));
    double second = 1. / (1. + exp(zetaOver2T * (x_u - R_u_k)));
    
    return first * second;
  }
};

// [[Rcpp::export]]
DataFrame hypermap(DataFrame net, double gma, double T, int k_speedup, double m_in, double L_in, double window, DataFrame theta){
  
  double zeta = 1.0; // zeta, to control the hyperbolic space curvature
  bool corrections = false; // to decide whether corrections steps should be performed (this is turned off to improve efficiency)
  
  // Read variables and compute L and m from kbar and min degree
  // Nodes are automatically sorted by degree
  Graph * G = new Graph();
  if(m_in != -1){ G->setM(m_in); } // if m has been specified by the user
  if(L_in != -1){ G->setL(L_in); } // if L has been specified by the user
  
  bool okRead = G->readGraph(net);
  if(!okRead){
    return NULL;
  }
  
  okRead = G->readLaBNEAngles(theta);
  if(!okRead){
    return NULL;
  }
  
  // Some variables
  const double N = G->getN();
  const double L = G->getL();
  const double m = G->getM();
  const double beta = 1. / (gma - 1);
  const double oneOverOneMinusBeta = 1. / (1. - beta);
  const double zetaOver2T = zeta / (2. * T);
  
  /*****************************************************************
   ******* Section of code that uses common neighbors begins *******
   ****************************************************************/
  
  int numCN = 0; // number of nodes which will be included in the CN computation
  
  // constant part of L_t (Eq. 3)
  const double cteL_t = 2.*L*(1.-beta) / (pow((1-pow(G->getN(),-(1-beta))),2) * (2*beta-1));
  
  // constant part or R_i (Eq. 2)
  const double cteR = 2./zeta * log(2.*T / (sin(T * M_PI)));
  
  // Compute R_t for each node
  vector<Node*> * nodes = G->getNodeList();
  
  for(int i = 0; i < N; ++i){
    
    double t = i+1;
    
    double I_t = oneOverOneMinusBeta * (1 - pow(t, -(1-beta))); // I_t
    
    double L_t;
    if(beta == 1)
      L_t = 2*L*(N-t)*log(t) / (t*pow(log(N),2));
    else if(beta == .5)
      L_t = L * ((1-pow(t,-0.5))/pow((1-pow(N,-0.5)),2)) * log(N/t);
    else
      L_t = cteL_t * (pow(N/t, (2*beta-1))-1) * (1-pow(t, -(1-beta)));
    
    double m_t = m + L_t;
    double r_t = (2./zeta) * log(t);
    double R_t = r_t - cteR - (2./zeta)*log(I_t/m_t);
    nodes->at(i)->setInitRadius(r_t, zeta);
    nodes->at(i)->setRadius(beta*r_t + (1-beta)*(2./zeta)*log(N), zeta);
    nodes->at(i)->setR(R_t);
    
    // Condition for including the node in the CN computation
    if((m_t >= t-1) || (t == 1))
      numCN++;
  }
  
  
  /***** REPLAY THE GROWTH *****/
  for(int i = 0; i < numCN; ++i){
    
    int t = i+1;
    
    Node *v = nodes->at(i);
    
    v->set_is_here();
    double labne_theta = v->getAngle();
    
    /* if(t == 1){
     // First node's angle is set to pi (it could be any other random value in [0, 2*pi]).
     v->setAngle(M_PI);
    }
     else{*/
    
    // Start calculating angles
    double step = min(1./t, 0.01);
    double maxLogL = -9999999999;
    
    //double theta_v = 0.00000001;
    double theta_v = labne_theta - (window/2);
    
    //while(theta_v <= 2*M_PI){
    while(theta_v <= labne_theta + (window/2)){
      double logL = 0.0;
      
      for(int j = 0; j < i; ++j){
        
        Node *u = nodes->at(j);
        double theta_u = u->getAngle();
        
        // Find the empirical # of common neighbors between v and u
        int empiricalCN = intersection(*(v->getAdjNodeList()), *(u->getAdjNodeList()));
        
        // Find the expected # of common neighbors between v and u
        double lambda = 0.0;
        double var = 0.0;
        
        for(int k = 0; k < N; ++k){
          if(k != i && k != j){
            
            Node *l = nodes->at(k);
            double r_l = l->getInitRadius();
            double r_v_k_1, r_v_k_2, rvk1_sinh, rvk2_sinh, rvk1_cosh, rvk2_cosh;
            double r_u_k_1, r_u_k_2, ruk1_sinh, ruk2_sinh, ruk1_cosh, ruk2_cosh;
            double R_u_k, R_v_k;
            if(r_l > v->getInitRadius()){ // l came after v
              r_v_k_2 = r_l;
              R_v_k = l->getR();
              r_v_k_1 = beta * v->getInitRadius() + (1-beta) * r_l;
              rvk1_sinh = sinh(zeta * r_v_k_1);
              rvk1_cosh = cosh(zeta * r_v_k_1);
              rvk2_sinh = l->getInitsinh();
              rvk2_cosh = l->getInitcosh();
            }
            else{ // v came after l
              r_v_k_2 = v->getInitRadius();
              R_v_k = v->getR();
              r_v_k_1 = beta * r_l + (1-beta) * v->getInitRadius();
              rvk1_sinh = sinh(zeta * r_v_k_1);
              rvk1_cosh = cosh(zeta * r_v_k_1);
              rvk2_sinh = v->getInitsinh();
              rvk2_cosh = v->getInitcosh();
            }
            
            if(r_l > u->getInitRadius()){ // l came after u
              r_u_k_2 = r_l;
              R_u_k = l->getR();
              r_u_k_1 = beta * u->getInitRadius() + (1-beta) * r_l;
              ruk1_sinh = sinh(zeta * r_u_k_1);
              ruk1_cosh = cosh(zeta * r_u_k_1);
              ruk2_sinh = l->getInitsinh();
              ruk2_cosh = l->getInitcosh();
            }
            else{ // u came after l
              r_u_k_2 = u->getInitRadius();
              R_u_k = u->getR();
              r_u_k_1 = beta * r_l + (1-beta) * u->getInitRadius();
              ruk1_sinh = sinh(zeta * r_u_k_1);
              ruk1_cosh = cosh(zeta * r_u_k_1);
              ruk2_sinh = u->getInitsinh();
              ruk2_cosh = u->getInitcosh();
            }
            
            // Numerical integration
            toIntegrate f = {zeta, zetaOver2T, theta_u, theta_v,
                             ruk1_sinh, ruk1_cosh, ruk2_sinh, ruk2_cosh,
                             rvk1_sinh, rvk1_cosh, rvk2_sinh, rvk2_cosh,
                             R_u_k, R_v_k};
            double prob = (1./(2*M_PI)) * stdfin::gauss_legendre_48(f, 0, 2*M_PI);
            
            //lambda is the total average (expected # of CN).
            lambda = lambda + prob;
            
            //The variance of each event is (1-prob)*prob and var is the total variance.
            var += (1-prob) * prob;
          }
        } //Endfor, for all possible neighbors k
        
        //Log of Normal distribution with mean lambda and variance var.
        double logProbCN = -pow(empiricalCN-lambda,2) / (2*var) - log(sqrt(2*M_PI*var));
        
        if(logProbCN > 0) {//if var is very small, since -log(sqrt(2*M_PI*var)) will become large.
          if((empiricalCN > lambda+1) || (empiricalCN < lambda-1)) //If empiricalCN far from lambda.
            logProbCN = -1000000000;
          else  //If empiricalCN close to lambda.
            logProbCN = 0;
        }
        logL += logProbCN;
      } //Endfor, for all nodes.
      
      if (logL > maxLogL) {
        v->setAngle(theta_v);
        maxLogL = logL;
      }
      
      theta_v += step;
      
    }//Endwhile, for all angles.
    
  }//Endfor, for all nodes.
  
  
  /***************************************************************
   ******* Section of code that uses common neighbors ends *******
   ***************************************************************/
  
  // Note: k_speedup should be smaller than the minimum degree for which we run correction steps,
  // otherwise we may alter the angular coordinates of existing nodes with approximate MLE
  // coordinates.
  // Setting k_speedup=0 disables the fast approximation.
  // Unless specified by the user, k_speedup=0.
  
  for(int i = numCN; i < N; ++i){
    
    int t = i+1;
    
    // R_t for all nodes is already calculated
    
    Node *u = nodes->at(i);
    
    u->set_is_here();
    u->setRadius(u->getInitRadius(), zeta);
    double labne_theta = u->getAngle();
    
    /*if(t == 1){
     // If first node, set its angle to pi (it can be any other random value in [0, 2*pi]).
     u->setAngle(M_PI);
     continue;
    }*/
    
    for(int j = 0; j < i; ++j){
      Node *v = nodes->at(j);
      v->setRadius(beta * v->getInitRadius() + (1-beta) * u->getInitRadius(), zeta);
    }
    
    vector<Node*> * nodes2compare;
    if(u->getK() < k_speedup)
      nodes2compare = u->getAdjNodeList(); // Consider only connections to my neighbors (that are already here).
    else
      nodes2compare = new vector<Node*>(nodes->begin(), nodes->begin() + i); // all older nodes
    
    // Start calculating angles
    double step = min(1./t, 0.01);
    double maxLogL = -9999999999;
    
    double theta = labne_theta - (window/2);
    
    while(theta <= labne_theta + (window/2)){
      
      double theta_t = theta;
      double logL = 0;
      
      for(int j = 0; j < nodes2compare->size(); ++j){
        
        Node *v = nodes2compare->at(j);
        
        // The check below is needed if we run the speedup heuristic as some of
        // the neighbors of v in the input graph may not have appeared yet.
        if (v->get_is_here() != 1) continue;
        
        double theta_s = v->getAngle();
        double dtheta = M_PI - fabs(M_PI - fabs(theta_t - theta_s));
        
        double x_st;
        if(dtheta == 0)
          x_st = fabs(v->getRadius() - u->getRadius());
        else
          x_st = (1./zeta) * acosh(u->getRcosh()*v->getRcosh() - u->getRsinh()*v->getRsinh()*cos(dtheta));
        
        double P_st = 1. / (1. + exp(zetaOver2T * (x_st - u->getR())));
        
        if(u->isNeighbor(v))
          logL += log(P_st);
        else
          logL += log(1-P_st);
        
      } // Endfor int j
      
      if(logL >= maxLogL){
        u->setAngle(theta_t);
        maxLogL = logL;
      }
      
      theta += step;
    } //Endwhile
    
    // Free memory
    if(u->getK() >= k_speedup)
      delete nodes2compare;
    
    // Only for k < k_speedup
    if(u->getK() < k_speedup){
      
      double C = 200.0;
      double theta_maxL = u->getAngle();
      double Delta = C * step;
      double theta_min = theta_maxL - Delta;
      double theta_max = theta_maxL + Delta;
      
      theta_min = max(theta_min, 0.0);
      theta_max = min(theta_max, 2*M_PI);
      
      maxLogL = -99999999999;
      
      double theta = theta_min;
      
      while(theta <= theta_max){
        
        double theta_t = theta;
        double logL = 0.0;
        
        for(int j = 0; j < i; ++j){
          
          Node *v = nodes->at(j);
          
          double theta_s = v->getAngle();
          double dtheta = M_PI - fabs(M_PI - fabs(theta_t - theta_s));
          
          double x_st;
          if(dtheta == 0)
            x_st = fabs(v->getRadius() - u->getRadius());
          else
            x_st = (1./zeta) * acosh(u->getRcosh()*v->getRcosh() - u->getRsinh()*v->getRsinh()*cos(dtheta));
          
          double P_st = 1. / (1. + exp(zetaOver2T * (x_st - u->getR())));
          
          if(u->isNeighbor(v))
            logL += log(P_st);
          else
            logL += log(1-P_st);
        }
        
        if(logL >= maxLogL){
          u->setAngle(theta_t);
          maxLogL = logL;
        }
        
        theta += step;
      } // Endwhile
    } // Endif
    
    /****************************
     ******* Corrections  *******
     ***************************/
    
    if(corrections){
      
      if(i+1 == N) continue;
      
      if (nodes->at(i+1)->getK() < u->getK() &&
          ((u->getK() == 60)
             || (u->getK() == 40)
             || (u->getK() == 20)
             || (u->getK() == 10))) {
             
             // We repeat each correction step a number of times, equal here to the average node degree
             // (fewer times is also beneficial.)
             for (int round = 1 ; round <= G->getKbar(); round++) {
               
               // consider updating the angular coordinates of all current nodes that were not inferred using common neighbors.
               vector<Node*> * nodes2compare = new vector<Node*>(nodes->begin()+numCN, nodes->begin()+i+1);
               
               // all current nodes.
               vector<Node*> * nodes2compare_all = new vector<Node*>(nodes->begin(), nodes->begin()+i+1);
               
               for(int j = 0; j < nodes2compare->size(); ++j){
                 
                 Node *v = nodes2compare->at(j);
                 
                 double step = min(1./t, 0.01);
                 double theta = 0.0;
                 double maxLogL = -99999999999;
                 
                 while(theta <= 2*M_PI){
                   
                   double theta_v = theta;
                   double logL = 0;
                   
                   for(int k = 0; k < nodes2compare_all->size(); ++k){
                     
                     Node *l = nodes2compare_all->at(k);
                     
                     if(v != l){
                       
                       double theta_l = l->getAngle();
                       double dtheta = M_PI - fabs(M_PI - fabs(theta_v - theta_l));
                       
                       double R_init;
                       double r_v, r_l;
                       
                       if(v->getInitRadius() > l->getInitRadius()){ // v came after l
                         r_v = v->getInitRadius();
                         r_l = beta * l->getInitRadius() + (1-beta) * r_v;
                         R_init = v->getR();
                       }
                       else{ // l came after v
                         r_l = l->getInitRadius();
                         r_v = beta * v->getInitRadius() + (1-beta) * r_l;
                         R_init = l->getR();
                       }
                       
                       double x_vl;
                       if(dtheta == 0)
                         x_vl = fabs(r_v - r_l);
                       else
                         x_vl=(1./zeta)*acosh((cosh(zeta*r_v)*cosh(zeta*r_l))-(sinh(zeta*r_v)*sinh(zeta*r_l)*cos(dtheta)));
                       
                       double P_vl = 1. / (1. + exp(zetaOver2T * (x_vl - R_init)));
                       
                       if(v->isNeighbor(l))
                         logL += log(P_vl);
                       else
                         logL += log(1-P_vl);
                     }
                   }
                   
                   if(logL >= maxLogL){
                     v->setAngle(theta_v);
                     maxLogL = logL;
                   }
                   
                   theta += step;
                   
                 } // Endwhile
                 
               }// Endfor pick next node
               
               // Free memory
               delete nodes2compare;
               delete nodes2compare_all;
               
             } // Endfor round
      } // Endfor correction step
    } // if(corrections)
  }// Endfor int i
  
  DataFrame coords = G->getCoords();
  
  // Free memory and return result
  delete G;
  return(coords);
}
