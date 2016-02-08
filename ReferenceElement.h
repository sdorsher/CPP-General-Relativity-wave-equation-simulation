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

  vector<double> refNodeLocationsVec;
  vector<double> refNodeWeightsVec;
  vector<double> vandermondeMatrixVec;
  vector<double> dVdrVec;
  vector<double> derivativeMatrixVec;
  vector<double> liftVec;
  vector<double> liftDim;
  vector<double> Ddim;
  
 public:
  const vector<double>& getD(); //get derivative matrix
  vector<double> getDdim();
  double getDelem(int, int); //get the derivative matrix at two indices
  const vector<double> getLift(); //get lift matrix
  vector<double>& getLiftDim();
  vector<double> getr(); //get node locations
  vector<double> getw(); //get weights
  int getOrder(); //get order
};

inline const vector<double>& ReferenceElement::getD()
{//Returns the derivative matrix.
  return derivativeMatrixVec;
}

inline vector<double> ReferenceElement::getDdim()
{//returns the dimensions of the derivative matrix.
  return Ddim;
}

inline double ReferenceElement::getDelem(int i, int j)
{//Returns the derivative matrix at i j
  return derivativeMatrixVec[i][j];
}

inline vector<double> ReferenceElement::getr()
{//Returns the reference node locations.
  return refNodeLocationsVec;
}

inline vector<double> ReferenceElement::getw()
{//Returns the reference node weights.
  return refNodeWeightsVec;
}

inline int ReferenceElement::getOrder()
{//Returns the element order.
  return order;
}

inline const vector<double>& ReferenceElement::getLift()
{//Returns the lift matrix.
  return liftVec;
}

inline vector<double> ReferenceElement::getLiftDim()
{//returns the dimensions of the lift matrix.
  return liftDim;
}

#endif
