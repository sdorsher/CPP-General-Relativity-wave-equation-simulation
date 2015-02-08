#include <iostream>
#include "TNT2.h"


using namespace TNT;

int main()
{
  Array2D<double> fours(4,4,4.0);
  Array1D<double> nines(3,9.0);

  std::cout << "two arrays" << std::endl;
  output2D(fours);
  output1D(nines);

  std::cout << "their sqrt" << std::endl;
  output2D(sqrt(fours));
  output1D(sqrt(nines));
  /// sqrt function works

  std::cout << "----------------------" << std::endl;

  std::cout << "a diagonal matrix" << std::endl;
  

  Array2D<double> diag(4,4,0.0);
  Array1D<double> vals(4);
  Array2D<double> vecs(4,4);

  diag[0][0]=4.;
  diag[1][1]=2.;
  diag[2][2]=3.;
  diag[3][3]=1.;

  output2D(diag);
  
  JAMA::Eigenvalue<double> eig(diag);
  eig.getRealEigenvalues(vals);
  eig.getV(vecs);
  std::cout << "its eigenvalues as found by JAMA" << std::endl;
  output1D(vals);
  std::cout << "its eigenvectors as found by Jama" << std::endl;
  output2D(vecs);
  //getRealEigenvalues and getV sort from least to greatest just like lapack


  /*  output=doubvec+2.0;

  output1D(output);
  output1D(doubvec);

  output+=doubvec;

  output1D(output);
 

  output=1.0+doubvec;

  output1D(output);

  output= 3.0*doubvec;
  output1D(output);

  output = doubvec * 4.0;
  output1D(output);

  output=doubvec/5.0;
  output1D(output);

  output=2.0/output;
  output1D(output);

  output = doubvec-3.0;
  output1D(output);
  
  output = 3.0-doubvec;
  output1D(output);
  //---------------- checks out up to here --------

  Array2D<double> paulix(2,2);
  Array2D<double> pauliz(2,2);
  Array2D<double> identity(2,2);

  paulix[0][1]=1.0;
  paulix[1][0]=1.0;
  pauliz[0][0]=1.0;
  pauliz[1][1]=-1.0;
  identity[0][0]=1.0;
  identity[1][1]=1.0;
  
  Array2D<double> out(2,2);

  out=identity+1.0;
  output2D(out);
  out=1.0+identity;
  output2D(out);
  
  out=1.0-identity;
  output2D(out);
  out=identity-1.0;
  output2D(out);
  
  out = 2.0*identity;
  output2D(out);
  out = identity*2.0;
  output2D(out);
  
  out=3.0/out;
  output2D(out);
  out=out/1.5;
  output2D(out);

  output=matmult(pauliz,doubvec);
  output1D(output);
  output=matmult(doubvec,pauliz);
  output1D(output);

  //checks out up to here -------------------

  */

}
