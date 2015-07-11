#include "CharacteristicFlux.h"


//This object DOES NOT CALCULATE THE CHARACTERISTIC FLUX. That is i
//misnomer. It converts the separated linearized differential equation
//to its characteristic form, which supplies the matrices needed
//to calculate the characteristic flux. This calculation is done
//in Grid::characteristicFlux 

// du/dt + A du/dx + Bu =0
// Atrimmed cuts the equation down to just those variables with nonzero 
// spatial derivatives 
// dw/dt + Lamb * dw/dx + ? =0, Lamb diagonal, characteristic form
// S matrix has eigenvectors of Atrimmed in columns
// Sinv= inverse of S
// Sinv*Atrimmed*S=Lamb, Lamb has eigenvalues of Atrimmed on diagonals

CharacteristicFlux::CharacteristicFlux(TNT::Array2D<double> Amatrix,
                                       TNT::Array2D<double> Atrim): 
  A(Amatrix.dim1(), Amatrix.dim2()),
  Atrimmed(Atrim.dim1(), Atrim.dim2()),
  Smatrix(0, 0), 
  Sinv(0, 0), 
  Lamb(0, 0), 
  one(0, 0)
{

  
  if(Amatrix.dim1() != Amatrix.dim2()) {
    throw invalid_argument("Amatrix not square in CharacteristicFlux");
  }

  if(Atrim.dim1() != Atrim.dim2()) {
    throw invalid_argument("Atrim not square in CharacteristicFlux");
  }

  A = Amatrix;
  Atrimmed = Atrim;
  Adimension = A.dim1();
  
  Ddimension = Atrim.dim1();

  Array2D<double> onetemp(Ddimension, Ddimension, 0.0);
  
  for(int i = 0; i < Ddimension; i++) {
    onetemp[i][i]=1.0;
  }
  
  one = onetemp.copy();
  
  Array2D<double> Smatrixtemp(Ddimension, Ddimension);
  Smatrix = Smatrixtemp.copy();

  JAMA::Eigenvalue<double> Seigen(Atrimmed);
  Seigen.getV(Smatrix);
  JAMA::LU<double> Sinverter(Smatrix);
  Sinv=Sinverter.solve(one);
  Lamb = matmult(Sinv, matmult(Atrimmed, Smatrix));

}

//this function does not work because it is not generally solveable. It works 
//for the wave equation but not in general
/*Array2D<double> CharacteristicFlux::trimA() 
{
  vector<int> non_zero_row_index;

  for(int i = 0; i < Adimension; i++) {//loop over rows
    bool rownonzeros = false;
    for( int j = 0; j < Adimension; j++) {//loop over columns 
      rownonzeros = ((A[i][j] != 0) || rownonzeros);
    }
    if(rownonzeros){
      non_zero_row_index.push_back(i); //save index of row
    }
  }

  Ddimension=non_zero_row_index.size();
  Array2D<double> Atrim(Ddimension, Ddimension);
  for(int i = 0; i < Ddimension; i++) {
    for(int j = 0; j < Ddimension; j++){
      Atrim[i][j] = A[nonzeros[i]][nonzeros[j]]; 
      //build Atrim out of nonzero rows
    }
  }
  return Atrim;
  }*/

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
