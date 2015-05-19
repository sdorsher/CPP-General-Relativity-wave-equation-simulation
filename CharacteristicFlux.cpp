#include "CharacteristicFlux.h"

CharacteristicFlux::CharacteristicFlux(TNT::Array2D<double> Amatrix): A(Amatrix.dim1(),Amatrix.dim2()),Smatrix(Amatrix.dim1(),Amatrix.dim1()), Sinv(Amatrix.dim1(),Amatrix.dim1()), Lamb(Amatrix.dim1(),Amatrix.dim1()), one(Amatrix.dim1(),Amatrix.dim1(),0.0)
{
  bool rowallzeros=true;
  int toptrimcount=0;
  int bottomtrimcount=0;

  A=Amatrix;

  if(Amatrix.dim1()!=Amatrix.dim2())
    {
      throw invalid_argument("Amatrix not square in characteristicFlux");
    }

  for(int i=0; i<one.dim1(); i++)
    {
      one[i][i]=1.0;
    }


  JAMA::Eigenvalue<double> Seigen(Amatrix);
  Seigen.getV(Smatrix);
  JAMA::LU<double> Sinverter(Smatrix);
  Sinv=Sinverter.solve(one);
  Lamb=matmult(Sinv,matmult(Amatrix,Smatrix));


}

Array2D<double> CharacteristicFlux::getS()
{
  return Smatrix;
}

Array2D<double> CharacteristicFlux::getA()
{
  return A;
}
Array2D<double> CharacteristicFlux::getSinv()
{
  return Sinv;
}
Array2D<double> CharacteristicFlux::getLambda()
{
  return Lamb;
}
