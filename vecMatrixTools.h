#ifndef VEC_MATRIX_TOOLS_H
#define VEC_MATRIX_TOOLS_H
#include <vector>

using namespace std;

vector<double> vecsum(vector <double> &A, vector<double> &B, int dimAB1, int dimAB2){
  vector<double> C(dimAB1*dimAB2);
  
  for(int i=0; i<dimAB1; i++){
    for(int j=0; j<dimAB2; j++){
      C[j*dimAB1+i]=A[j*dimAB1+i]+B[j*dimAB1+i];
    }
  }
  return C;
}

vector<double> vecdiff(vector <double> &A, vector<double> &B, int dimAB1, int dimAB2){
  vector<double> C(dimAB1*dimAB2);
  
  for(int i=0; i<dimAB1; i++){
    for(int j=0; j<dimAB2; j++){
      C[j*dimAB1+i]=A[j*dimAB1+i]-B[j*dimAB1+i];
    }
  }
  return C;
}

vector<double> scalarmult(double s, vector<double> &A, int dimA1, int dimA2){
  vector<double> C(dimA1*dimA2);
  
  for(int i=0; i<dimA1; i++){
    for(int j=0; j<dimA2; j++){
      C[j*dimA1+i]=A[j*dimA1+i]*s;
    }
  }
  return C;
}

vector<double> matmul(vector<double> &A, vector<double> &B, int dimA1,
		       int dimA2B1, int dimB2){
  vector<double> C(dimA1*dimB2);
  
  for(int i=0; i<dimA1; i++){
    for(int j=0; j<dimB2; j++){
      double sum = 0.0;
      for(int k=0; k<dimA2B1;k++){
	sum+=A[i*dimA2B1+k]*B[k*dimB2+j];
      }
      C[i*dimB2+j]=sum;
    }
  }
  return C;

}
#endif
