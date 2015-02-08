#include "tnt_array1d.h"
#include "tnt_array2d."
#include "tnt_array1d_utils.h"
#include "tnt_array2d_utils.h"

class ReferenceElement
{
 private:
  int order; //order of element
  Array1D refNodeLocations; // node locations scaled to r
  Array2D vandermondeMatrix; 
  Array2D derivativeMatrix;
  Array2D dVdr;
 public:
  void initialize(int N); //static initializer outside of constructor
 private:
  void jacobiPolynomails();
  void jglQuadraturePoints();
  void computeVandermonde();
  void computedVdr();
  void computeDr();
 public:
  Array2D getDr(); //get derivative matrix
};
  
