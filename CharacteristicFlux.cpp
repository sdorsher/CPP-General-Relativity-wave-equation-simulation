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
// S matrix has eigevectors of Atrimmed in columns
// Sinv= inverse of S
// Sinv*Atrimmed*S=Lamb, Lamb has eigenvalues of Atrimmed on diagonals

//We can leave TNT in here because it is only done at the beginning

CharacteristicFlux::CharacteristicFlux(vector<double> Amatrix,
                                       vector<double> Atrim)
  //  A(Amatrix.dim1(), Amatrix.dim2()),
  //  Atrimmed(Atrim.dim1(), Atrim.dim2()),
  //  Smatrix(0, 0), 
  //  Sinv(0, 0), 
  //  Lamb(0, 0), 
  //  one(0, 0)
{

  
  //  if(Amatrix.dim1() != Amatrix.dim2()) {
  //  throw invalid_argument("Amatrix not square in CharacteristicFlux");
  // }

  //if(Atrim.dim1() != Atrim.dim2()) {
  //  throw invalid_argument("Atrim not square in CharacteristicFlux");
  // }

  
  Array2D<double> A(params.grid.Adim,params.grid.Adim);
  Array2D<double> Atrimmed(params.grid.Ddim,params.grid.Ddim);
  Array2D<double> Smatrix(0, 0); 
  Array2D<double> Sinv(0, 0); 
  Array2D<double> Lamb(0, 0); 
  Array2D<double> one(0, 0,1.0);
  A = vectorToArray2D(Amatrix,params.grid.Adim,params.grid.Adim);
  Atrimmed = vectorToArray2D(Atrim,params.grid.Ddim,params.grid.Ddim);
  Adimension = params.grid.Adim;
  Ddimension = params.grid.Ddim;

  Array2D<double> onetemp(Ddimension, Ddimension, 0.0);
  
  for(int i = 0; i < Ddimension; i++) {
    onetemp[i][i]=1.0;
  }
  
  one = onetemp.copy();
  
  Array2D<double> Smatrixtemp(Ddimension, Ddimension);
  Smatrix = Smatrixtemp.copy();
  JAMA::Eigenvalue<double> Seigen(Atrimmed);
  Seigen.getV(Smatrix); // get Eigenvector matrix and store it in Smatrix
  JAMA::LU<double> Sinverter(Smatrix); 
  Sinv=Sinverter.solve(one); //Invert the Eigenvector matrix to get Sinv
  Lamb = matmult(Sinv, matmult(Atrimmed, Smatrix));
  //Find the characteristic matrix for the system of equations, Lamb


  SmatrixV=Array2DtoVector(Smatrix);
  SinvV=Array2DtoVector(Sinv);
  LambV=Array2DtoVector(Lamb);
  AtrimmedV=Array2DtoVector(Atrimmed);
  AV=Array2DtoVector(A);

}

vector<double> CharacteristicFlux::getS(){
  return SmatrixV;
}

vector<double> CharacteristicFlux::getA(){
  return AV;
}

vector<double> CharacteristicFlux::getSinv(){
  return SinvV;
}

vector<double> CharacteristicFlux::getLambda(){
  return LambV;
}

vector<double> CharacteristicFlux::getAtrimmed(){
  return AtrimmedV;
}
  
int CharacteristicFlux::getAdim(){
  return Adimension;
}

int CharacteristicFlux::getDdim(){
  return Ddimension;
}
