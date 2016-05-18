#ifndef VEC_MATRIX_TOOLS_H
#define VEC_MATRIX_TOOLS_H
#include <vector>
#include <iostream>

using namespace std;

template <class T>
void insert_1D_into_2D_vec(vector<T> &A, vector<T> &B, int dimA1,int dimA2, int index, bool iscolumn){
  if(iscolumn){
    if(index>=dimA2) cout << "index out of range in insert_1D_into_2D_vec" << endl;
    for(int j=0; j<dimA2; j++){
      A[index*dimA2+j]=B[j];
    }
  }else{
    if(index>=dimA1) cout << "index out of range in insert_1D_into_2D_vec" << endl;
    for(int i=0; i<dimA1; i++){
      A[i*dimA2+index]=B[i];
    }
  }
}

template <class T>
vector<T> vecsum(vector <T> &A, vector<T> &B){
  vector<T> C(A.size());

  if(A.size()!=B.size()) cout << "vecsum dimension mismatch\n";
  
  for(int i=0; i<A.size(); i++){
    C[i]=A[i]+B[i];
  }
  return C;
}

template <class T>
vector<T> vecdiff(vector <T> &A, vector<T> &B){
  vector<T> C(A.size());

  if(A.size()!=B.size()) cout << "vecdiff dimension mismatch\n";
  for(int i=0; i<A.size(); i++){
    C[i]=A[i]-B[i];
  }
    return C;
}

template <class T>
vector<T> scalarmult(T s, vector<T> &A){
  vector<T> C(A.size());
  
  for(int i=0; i<A.size(); i++){
      C[i]=A[i]*s;
    }
  
  return C;
}

template <class T>
vector<complex<T>> scalarmult(complex<T> s, vector<T> &A){
  vector<complex<T>> C(A.size());
  
  for(int i=0; i<A.size(); i++){
      C[i]=A[i]*s;
    }
  
  return C;
}

template <class T>
vector<complex<T>> scalarmult(T s, vector<complex<T>> &A){
  vector<complex<T>> C(A.size());
  
  for(int i=0; i<A.size(); i++){
      C[i]=A[i]*s;
    }
  
  return C;
}


template <class T>
vector<T> matmul(vector<T> &A, vector<T> &B, int dimA1,
		       int dimA2B1, int dimB2){
  vector<T> C(dimA1*dimB2);
  
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

template <class T>
vector<complex<T>> matmul(vector<T> &A, vector<complex<T>> &B, int dimA1,
		       int dimA2B1, int dimB2){
  vector<complex<T>> C(dimA1*dimB2);
  
  for(int i=0; i<dimA1; i++){
    for(int j=0; j<dimB2; j++){
      complex<double> sum = 0.0;
      for(int k=0; k<dimA2B1;k++){
	sum+=A[i*dimA2B1+k]*B[k*dimB2+j];
      }
      C[i*dimB2+j]=sum;
    }
  }
  return C;

}
#endif
