#include "../tnt/tnt.h"
#include "../jama125/jama_lu.h"
#include "../jama125/jama_eig.h"
using namespace std;

void output2D(TNT::Array2D<double> matr)
{
  for(int i=0;i<matr.dim1();i++){
    cout << "( ";
    for(int j=0;j<matr.dim2();j++){
      cout << matr[i][j] <<" "; 
    }
    cout << ")\n";
  }
  cout << endl << endl;
}

void output1D(TNT::Array1D<double> matr)
{
  for(int i=0;i<matr.dim();i++){
    cout << "( ";
    cout << matr[i] <<" "; 
    cout << ")\n";
  }
  cout << endl<< endl;
}


int main()
{
  TNT::Array2D<double> vec(3,1);
  TNT::Array2D<double> matr1(2,3);
  TNT::Array2D<double> matr2(2,2);
  TNT::Array2D<double> matr3(2,3);
  TNT::Array2D<double> matr4(2,2);
  TNT::Array2D<double> ident(2,2);

  vec[0][0]=1;
  vec[1][0]=2;
  vec[2][0]=3;

  matr1[0][0]=1;
  matr1[0][1]=1;
  matr1[0][2]=2;
  matr1[1][0]=2;
  matr1[1][1]=0;
  matr1[1][2]=3;
 
  matr2[0][0]=0;
  matr2[0][1]=1;
  matr2[1][0]=1;
  matr2[1][1]=0;

  matr3[0][0]=-1;
  matr3[0][1]=-1;
  matr3[0][2]=-1;
  matr3[1][0]=1;
  matr3[1][1]=1;
  matr3[1][2]=1;

  matr4[0][0]=1;
  matr4[0][1]=2;
  matr4[0][2]=-1;
  matr4[0][1]=3;

  ident[0][0]=1;
  ident[1][0]=0;
  ident[0][1]=0;
  ident[1][1]=1;

  cout << matr1.dim1() << "\t" << matr1.dim2() <<"\t" << vec.dim1() << endl;

  cout << "matr1" << endl;
  output2D(matr1);
  cout << "matr2" <<endl;
  output2D(matr2);
  cout << "matr3" <<endl;
  output2D(matr3);
  cout << "vec" << endl;
  output2D(vec);
  
  cout << "==========" << endl;

  //must use matmult instead of times
  TNT::Array2D<double> prodmatr=TNT::matmult(matr2,matr1);
  TNT::Array2D<double> summatr = matr1+matr3;
  //no exponentiation. not sure that is what ** means anyway
  //TNT::Array2D<double> expmatr = matr2**3;
  TNT::Array2D<double> prodvec = TNT::matmult(matr1,vec);
  //no matrix multiplication between 1D and 2D array
  //no scalar multiplication?
  //TNT::Array2D<double> scalarmultmatr = 2.0*matr2;
  //no insertion of 1D array into 2D array, or smaller 2D array into 2D array
  

  cout << "matr2*matr1" << endl;
  output2D(prodmatr);
  
  cout << "matr1+matr3" << endl;
  output2D(summatr);
  
  //cout<< "matr2**3" << endl;
  //output2D(expmatr);
  
  cout<< "matr1*vec" << endl;
  output2D(prodvec);
  
  
  //  cout << "2.0*matr2" <<endl;
  // output2D(scalarmultmatr);

  cout << "================" <<endl;

  cout << "original vector" <<endl;
  output2D(vec);
  TNT::Array2D<double> movevec = vec;
  TNT::Array2D<double> copyvec = vec.copy();
  vec[2][0]=2;
  cout << "original modified vector" <<endl;
  output2D(vec); 
  cout << "moved also changed" <<endl;
  output2D(movevec);
  cout << "copied not changed" <<endl;
  output2D(copyvec);

  cout << "==============" <<endl;
  
  cout << "matr4" <<endl;
  
  output2D(matr4);

  cout << "test eigenvalues matr4" << endl;

  JAMA::Eigenvalue<double> eig4(matr4);
  TNT::Array1D<double> realeig4;
  TNT::Array1D<double> imageig4;
  eig4.getRealEigenvalues(realeig4);
  eig4.getImagEigenvalues(imageig4);
  cout << "real part of eigenvalues" << endl;
  output1D(realeig4);
  cout << "imaginary part of eigenvalues " << endl;
  output1D(imageig4);
  

  cout << "===============" <<endl;

  cout << "test inverting matr4" << endl;

  JAMA::LU<double> matrtoinv4(matr4);
  TNT::Array2D<double> invmatr4 = matrtoinv4.solve(ident);

  output2D(invmatr4);

}

