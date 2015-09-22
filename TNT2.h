
#ifndef TNT2_H
#define TNT2_H

#include "tnt/tnt.h"
#include <typeinfo>
#include <iostream>
#include <stdexcept>
#include "tnt_array2D_extn.h"
#include "tnt_array1D_extn.h"
#include "jama/jama_eig.h"
#include "jama/jama_lu.h" 
#include <complex>


using namespace std;
namespace TNT
{

  //multiplies a complex Array2D by a complex Array2D
  template <class T>
    Array2D<complex<T>> matmult(const Array2D<complex<T>> &A, const Array2D<T> &B)
    {
      if (A.dim2() != B.dim1())
        return Array2D<complex<T>>();
      
      int M = A.dim1();
      int N = A.dim2();
      int K = B.dim2();
      
      Array2D<complex<T>> C(M,K);
    
      for (int i=0; i<M; i++)
        for (int j=0; j<K; j++){
	  complex<T> sum = 0;
	  
	  for (int k=0; k<N; k++)
	    sum += A[i][k] * B [k][j];
	    
	  C[i][j] = sum;
	}
      
      return C;
    }

  //Multiplies a real Array2D by a complex Array2D to the right
  template <class T>
    Array2D<complex<T>> matmult(const Array2D<T> &A, const Array2D<complex<T>> &B)
    {
      if (A.dim2() != B.dim1())
        return Array2D<complex<T>>();
      
      int M = A.dim1();
      int N = A.dim2();
      int K = B.dim2();
    
      Array2D<complex<T>> C(M,K);
      
      for (int i=0; i<M; i++)
        for (int j=0; j<K; j++){
	  complex<T> sum = 0;
	  
	  for (int k=0; k<N; k++)
	    sum += A[i][k] * B [k][j];
	  
	  C[i][j] = sum;
        }
      
      return C;
    }

  //Multiplies a 2D matrix by a 1D column vector to the right.
  template <typename T>
    TNT::Array1D<T> matmult(const TNT::Array2D<T>& A, const TNT::Array1D<T>& B)
    {
      if(A.dim2()!=B.dim1()){
        throw std::invalid_argument("Mismatched matrix dimensions.");
      }
    
      TNT::Array1D<T> ans(A.dim1());
    
      for(int i=0; i<A.dim1(); i++){
        ans[i]=0;
        for(int j=0; j<A.dim2(); j++){
          ans[i]+=A[i][j]*B[j];
        }
      }
      return ans;
    }

  //Multiplies a 1D row vector by a 2D matrix
  template <typename T>
    TNT::Array1D<T> matmult(const TNT::Array1D<T>& A, const TNT::Array2D<T>& B)
    {//really using the transpose of A, produces a column vector when it 
      //should produce a row vector
      if(A.dim1()!=B.dim1()){
        throw std::invalid_argument("Mismatched matrix dimensions.");
      }
    
      TNT::Array1D<T> ans(B.dim2());
    
      for(int j=0; j<B.dim2(); j++){
        ans[j]=0;
        for(int i=0; i<B.dim1(); i++){
          ans[j]+=A[i]*B[i][j];
        }
      }
      return ans;
    }
  
  //Inserts a 1D array into a 1D array at index origin
  template<typename T>
    void insert_1D(TNT::Array1D<T>& original, 
                   const TNT::Array1D<T>& inserted, const int origin)
    {
      if(origin+inserted.dim()>original.dim()){
        throw std::out_of_range("Inserted array dimensions plus origin exceed \n original array dimensions.");
      }
      for(int j=origin; j<origin+inserted.dim(); j++){
        original[j]=inserted[j-origin];
      }
    }
  
  //Inserts a 2D array into a 2D array at indices origin1,origin2
  template <typename T>
    void insert_2D(TNT::Array2D<T>& original, const TNT::Array2D<T>& inserted, 
		   const int origin1, const int origin2)
    //origin gives the index of the upper left most element of the 
    //inserted array
    {
      if((origin1+inserted.dim1()>original.dim1())||
	 (origin2+inserted.dim2()>original.dim2())){
        throw std::out_of_range("Inserted array dimensions plus origin exceed \n original array dimensions.");
      }
      for(int i=origin1; i<origin1+inserted.dim1(); i++){
        for(int j=origin2; j<origin2+inserted.dim2(); j++){
          original[i][j]=inserted[i-origin1][j-origin2];
        }
      }
    }
  
  //Inserts a 1D array into a 2D array at index along either a column or a 
  //row depending on whether iscolumn is true or false. 
  template <typename T>
    void insert_1D_into_2D(TNT::Array2D<T>& original, 
			   const TNT::Array1D<T>& inserted, 
			   const int index, bool iscolumn)
    //index gives the row or column into which the array should 
    //be inserted, iscolumn gives whether the 1D array is a 
    //column or row vector
    {
      if(iscolumn){
	if((index>=original.dim2())
	   ||(index<0)
	   ||(inserted.dim()!=original.dim1())){
          throw std::out_of_range("Inserted array dimension not equal to external array dimension  or index out of range.");
        }
	
	for(int i=0; i<inserted.dim(); i++){
          original[i][index]=inserted[i];
        }
      } else{//is row
	if((index>=original.dim1())
	   ||(index<0)
	   ||(inserted.dim()!=original.dim2()))
	  {
	    throw std::out_of_range("Inserted array dimension not equal to external array dimension or index out of range.");
	  }
        for(int i=0; i<inserted.dim(); i++){
          original[index][i]=inserted[i];
        }
      }    
    }
 
  //prints an Array1D
  template <typename T>
    void output1D(const TNT::Array1D<T>& matr)
    {
      for(int i=0;i<matr.dim();i++){
        cout << "( ";
        cout << matr[i] <<" "; 
        cout << ")\n";
      }
      cout << endl<< endl;
    }
  
  //prints an Array2D
  template <typename T>
    void output2D(const TNT::Array2D<T>& matr)
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

  //prints an Array1D
  template <typename T>
    void output1Dcomplex(const TNT::Array1D<complex<T>>& matr)
    {
      for(int i=0;i<matr.dim();i++){
        cout << "( ";
	  cout << matr[i].real() << "+" << matr[i].imag()<< "i" <<" "; 
        cout << ")\n";
      }
      cout << endl<< endl;
    }
  
  //prints an Array2D
  template <typename T>
    void output2Dcomplex(const TNT::Array2D<complex<T>>& matr)
    {
      for(int i=0;i<matr.dim1();i++){
        cout << "( ";
        for(int j=0;j<matr.dim2();j++){
	  cout << matr[i][j].real() << "+" << matr[i][j].imag()<< "i" <<" "; 
	}
        cout << ")\n";
      }
      cout << endl << endl;
    }

  //Multiplies a 2D matrix by a 1D column vector to the right.
  template <typename T>
    TNT::Array1D<complex<T>> matmult(const TNT::Array2D<T>& A, const TNT::Array1D<complex<T>>& B)
    {
      if(A.dim2()!=B.dim1()){
        throw std::invalid_argument("Mismatched matrix dimensions.");
      }
    
      TNT::Array1D<complex<T>> ans(A.dim1());
    
      for(int i=0; i<A.dim1(); i++){
        ans[i]=0;
        for(int j=0; j<A.dim2(); j++){
          ans[i]+=A[i][j]*B[j];
        }
      }
      return ans;
    }

  //Multiplies a 1D row vector by a 2D matrix
  template <typename T>
    TNT::Array1D<complex<T>> matmult(const TNT::Array1D<complex<T>>& A, const TNT::Array2D<T>& B)
    {//really using the transpose of A, produces a column vector when it 
      //should produce a row vector
      if(A.dim1()!=B.dim1()){
        throw std::invalid_argument("Mismatched matrix dimensions.");
      }
    
      TNT::Array1D<complex<T>> ans(B.dim2());
    
      for(int j=0; j<B.dim2(); j++){
        ans[j]=0;
        for(int i=0; i<B.dim1(); i++){
          ans[j]+=A[i]*B[i][j];
        }
      }
      return ans;
    }
  
  
};
  
#endif 
