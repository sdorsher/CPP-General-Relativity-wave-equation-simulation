#ifndef CHARACTERISTICFLUX_H
#define CHARACTERISTICFLUX_H
#include "TNT2.h"

using namespace TNT;
using namespace JAMA;

//du/dt + A du/dx +Bu=0

//dw/dt + Lamb dw/dx + ?=0 characteristic form

class CharacteristicFlux
{
 private:
  Array2D<double> A; //see above
  Array2D<double> Atrimmed; //part of A that has nonzero derivative dependence
  Array2D<double> Smatrix; //matrix with eigenvectors of A in columns
  Array2D<double> Sinv; //inverse of S matrix
  Array2D<double> Lamb; //eigenvalue of A matrix (characteristic form)
  Array2D<double> one; //identity matrix
  //  Array2D<double> trimA(); 
  int Adimension; //dimension of A
  int Ddimension; //dimension of non-zero derivative variables

 public:
  CharacteristicFlux(Array2D<double> Amatrix, Array2D<double> Atrimmed);
  Array2D<double> getA();
  Array2D<double> getS();
  Array2D<double> getSinv();
  Array2D<double> getLambda();
  Array2D<double> getAtrimmed();
  int getAdim();//untrimmed A dimension
  int getDdim();//derivative dimension (Atrimmed.dim())
};

#endif
