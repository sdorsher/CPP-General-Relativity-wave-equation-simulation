#include "CharacteristicFlux.h"

CharacteristicFlux::CharacteristicFlux(TNT::Array2D<double> Amatrix): 
  A(Amatrix.dim1(),Amatrix.dim2()),
  Atrimmed(0,0),
  Smatrix(0,0), 
  Sinv(0,0), 
  Lamb(0,0), 
  one(0,0)
{

  
  if(Amatrix.dim1()!=Amatrix.dim2())
    {
      throw invalid_argument("Amatrix not square in CharacteristicFlux");
    }

  A=Amatrix;
  Adimension=A.dim1();

  Atrimmed=trimA();
  
  Array2D<double> onetemp(Ddimension,Ddimension,0.0);
  
  for(int i=0; i<Ddimension; i++)
    {
      onetemp[i][i]=1.0;
    }
  
  one=onetemp.copy();

  Array2D<double> Smatrixtemp(Ddimension,Ddimension);
  Smatrix=Smatrixtemp.copy();
  
  JAMA::Eigenvalue<double> Seigen(Atrimmed);
  Seigen.getV(Smatrix);
  JAMA::LU<double> Sinverter(Smatrix);
  Sinv=Sinverter.solve(one);
  Lamb=matmult(Sinv,matmult(Atrimmed,Smatrix));


}

Array2D<double> CharacteristicFlux::trimA()
{
  vector<int> nonzeros;
   for(int i=0; i<Adimension; i++)
    {//loop over rows
      bool rownonzeros=false;
     for( int j=0; j<Adimension; j++)
        {//loop over columns
          rownonzeros=((A[i][j]!=0)||rownonzeros);
        }
      if(rownonzeros)
        {
          nonzeros.push_back(i);
        }
    }

  Ddimension=nonzeros.size();
  Array2D<double> Atrim(Ddimension,Ddimension);
  for(int i=0; i<Ddimension; i++)
    {
      for(int j=0; j<Ddimension; j++)
        {
          Atrim[i][j]=A[nonzeros[i]][nonzeros[j]];
        }
    }

  return Atrim;

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
Array2D<double> CharacteristicFlux::getAtrimmed()
{
  return Atrimmed;
}

  
int CharacteristicFlux::getAdim()
{
  return Adimension;
}

int CharacteristicFlux::getDdim()
{
  return Ddimension;
}
