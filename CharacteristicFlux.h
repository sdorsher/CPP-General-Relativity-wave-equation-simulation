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
  vector<double> AV; //see above
  vector<double> AtrimmedV; //part of A that has nonzero derivative dependence
  vector<double> SmatrixV; //matrix with eigenvectors of A in columns
  vector<double> SinvV; //inverse of S matrix
  vector<double> LambV; //eigenvalue of A matrix (characteristic form)
  //  Array2D<double> trimA(); 
  int Adimension; //dimension of A
  int Ddimension; //dimension of non-zero derivative variables

 public:
  CharacteristicFlux(Array2D<double> Amatrix, Array2D<double> Atrimmed);
  vector<double> getA();
  vector<double> getS();
  vector<double> getSinv();
  vector<double> getLambda();
  vector<double> getAtrimmed();
  int getAdim();//untrimmed A dimension
  int getDdim();//derivative dimension (Atrimmed.dim())
};

#endif
