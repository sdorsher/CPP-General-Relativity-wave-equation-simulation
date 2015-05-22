#include "TNT2.h"

using namespace TNT;
using namespace JAMA;

class CharacteristicFlux
{
 private:
  Array2D<double> A;
  Array2D<double> Atrimmed;
  Array2D<double> Smatrix;
  Array2D<double> Sinv;
  Array2D<double> Lamb;
  Array2D<double> one;
  Array2D<double> trimA();
  int Adimension;
  int Ddimension;

 public:
  CharacteristicFlux(Array2D<double> Amatrix);
  Array2D<double> getA();
  Array2D<double> getS();
  Array2D<double> getSinv();
  Array2D<double> getLambda();
  Array2D<double> getAtrimmed();
  int getAdim();//untrimmed A dimension
  int getDdim();//derivative dimension
};
