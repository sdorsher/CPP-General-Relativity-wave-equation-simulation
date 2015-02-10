#include <cmath>
#include "tnt.h"
#include "TNT2.h"

using namespace TNT;

class ReferenceElement
{

 public:
 ReferenceElement(int N);

  void jacobiGQ(TNT::Array1D<double>& x, double alpha, double beta, int n, 
		TNT::Array1D<double>& w);
  Array1D<double> jacobiGL(double alpha, double beta, double n);
  Array1D<double> jacobiP(const TNT::Array1D<double>& x, double alpha, 
	       double beta, int N);
  void vandermonde1D();

public:
  int order; //order of element
  Array1D<double> refNodeLocations; // node locations scaled to r
  Array2D<double> vandermondeMatrix; 
  //Array2D dVdr;
  //Array2D derivativeMatrix;

//Evaluate Jacobi Polynomial of type (alpha,beta) at points x for order N
  // void jglQuadraturePoints();
   //void computedVdr();
  //void computeDr();
  // public:
  //Array2D getDr(); //get derivative matrix
};
  
