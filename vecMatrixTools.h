#ifndef VEC_MATRIX_TOOLS_H
#define VEC_MATRIX_TOOLS_H
#include <vector>
#include <iostream>

using namespace std;

void insert_1D_into_2D_vec(vector<double> &A, vector<double> &B, int dimA1,int dimA2, int index, bool iscolumn){
  if(iscolumn){
    if(index>=dimA2) cout << "index out of range in insert_1D_into_2D_vec" << endl;
    for(int j=0; j<dimA2; j++){
      A[j*dimA1+index]=B[j];
    }
  }else{
    if(index>=dimA1) cout << "index out of range in insert_1D_into_2D_vec" << endl;
    for(int i=0; i<dimA1; i++){
      A[index*dimA1+i]=B[i];
    }
  }
}


vector<double> vecsum(vector <double> &A, vector<double> &B){
  vector<double> C(A.size());

  if(A.size()!=B.size()) cout << "vecsum dimension mismatch\n";
  
  for(int i=0; i<A.size(); i++){
    C[i]=A[i]+B[i];
  }
  return C;
}

vector<double> vecdiff(vector <double> &A, vector<double> &B){
  vector<double> C(A.size());

  if(A.size()!=B.size()) cout << "vecdiff dimension mismatch\n";
  for(int i=0; i<A.size(); i++){
    C[i]=A[i]-B[i];
  }
    return C;
}

vector<double> scalarmult(double s, vector<double> &A){
  vector<double> C(A.size());
  
  for(int i=0; i<A.size(); i++){
      C[i]=A[i]*s;
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
