#ifndef REFERENCE_ELEMENT_H
#define REFERENCE_ELEMENT_H
#include <cmath>
#include "tnt/tnt.h"
#include "TNT2.h"
#include "globals.h"

using namespace TNT;

class ReferenceElement
{

 public:
 ReferenceElement(int N);
 

 private:
  void jacobiGQ(TNT::Array1D<double>& x, double alpha, double beta, int n, 
		TNT::Array1D<double>& w);
  Array1D<double> jacobiGL(double alpha, double beta, double n);
  Array1D<double> gaussWeights(double alpha, double beta, int n);
  Array1D<double> jacobiP(const TNT::Array1D<double>& x, double alpha, 
	       double beta, int N);
  void vandermonde1D();
  Array1D<double> gradJacobiP(double alpha, double beta,int N);//evaluated at nodes
  void gradVandermonde1D(); //evaluated at nodes for order of element
  void Dmatrix1D();
  void lift1D();

private:
  int order; //order of element
  Array1D<double> refNodeLocations; // node locations scaled to r
  Array1D<double> refNodeWeights;
  Array2D<double> vandermondeMatrix; 
  Array2D<double> dVdr;
  Array2D<double> derivativeMatrix;
  Array2D<double> lift;

 public:
  Array2D<double> getD(); //get derivative matrix
  Array1D<double> getr(); //get node locations
  int getOrder(); //get order
  Array2D<double> getLift(); //get lift matrix
  Array1D<double> getw(); //get weights
  
};

//extern ReferenceElement refelem;
  
#endif
