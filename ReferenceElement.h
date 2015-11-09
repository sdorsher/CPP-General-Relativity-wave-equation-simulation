#ifndef REFERENCE_ELEMENT_H
#define REFERENCE_ELEMENT_H
#include <cmath>
#include "tnt/tnt.h"
#include "TNT2.h"
#include "globals.h"
#include <iomanip>

using namespace TNT;

class ReferenceElement
{

  //see .cpp file for explanation of functions
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
  void Dmatrix1D(); //calculate derivative matrix
  void lift1D(); //calculate lift matrix to be used in computation of flux

private:
  int order; //order of element
  Array1D<double> refNodeLocations; // node locations scaled to r
  Array1D<double> refNodeWeights;
  Array2D<double> vandermondeMatrix; 
  Array2D<double> dVdr;
  Array2D<double> derivativeMatrix;
  Array2D<double> lift; //used in numerical flux calculation, 
                        //scaled by jacobian

 public:
  Array2D<double> getD(); //get derivative matrix
  double getDelem(int, int); //get the derivative matrix at two indices
  Array2D<double> getLift(); //get lift matrix
  Array1D<double> getr(); //get node locations
  Array1D<double> getw(); //get weights
  int getOrder(); //get order
};

inline Array2D<double> ReferenceElement::getD()
{//Returns the derivative matrix.
  return derivativeMatrix;
}


inline double ReferenceElement::getDelem(int i, int j)
{//Returns the derivative matrix at i j
  return derivativeMatrix[i][j];
}

inline Array1D<double> ReferenceElement::getr()
{//Returns the reference node locations.
  return refNodeLocations;
}

inline Array1D<double> ReferenceElement::getw()
{//Returns the reference node weights.
  return refNodeWeights;
}

inline int ReferenceElement::getOrder()
{//Returns the element order.
  return order;
}

inline Array2D<double> ReferenceElement::getLift()
{//Returns the lift matrix.
  return lift;
}



#endif
