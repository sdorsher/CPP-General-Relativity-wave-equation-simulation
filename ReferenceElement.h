#include <cmath>
#include "../tnt/TNT.h"
#include "../tnt2/TNT2.h"


class ReferenceElement
{

 public:
  ReferenceElement(int N);

  void jacobiGQ(TNT::Array1D<double>& x, double alpha, double beta, int n, 
		TNT::Array1D<double>& w);
  Array1D<double> jacobiGL(double alpha, double beta, double n);
  void jacobiP(const TNT::Array1D<double>& x, double alpha, 
	       double beta, int N, TNT::Array1D<double>& polynom);
 private:
  int order; //order of element
  //Array1D refNodeLocations; // node locations scaled to r
  // Array2D vandermondeMatrix; 
  //Array2D dVdr;
  //Array2D derivativeMatrix;

//Evaluate Jacobi Polynomial of type (alpha,beta) at points x for order N
  // void jglQuadraturePoints();
  //void computeVandermonde();
  //void computedVdr();
  //void computeDr();
  // public:
  //Array2D getDr(); //get derivative matrix
};
  
