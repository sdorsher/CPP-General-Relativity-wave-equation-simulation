#ifndef TNT2_H
#define TNT2_H

#include "../tnt/tnt.h"
#include <typeinfo>
#include <iostream>
#include <stdexcept>
#include "tnt_array2d_extn.h"
#include "tnt_array1d_extn.h"
#include "../jama125/jama_eig.h"
#include "../jama125/jama_lu.h" //is this right? inverse, right?


using namespace std;
namespace TNT
{

  
  template <typename T>
    TNT::Array1D<T> matmult(const TNT::Array2D<T>& A, const TNT::Array1D<T>& B)
  {
    if(A.dim2()!=B.dim1())
      {
	throw invalid_argument("Mismatched matrix dimensions.");
      }
    TNT::Array1D<T> ans(A.dim1());
    
    for(int i=0; i<A.dim1(); i++)
      {
	ans[i]=0;
	for(int j=0; j<A.dim2(); j++)
	  {
	    ans[i]+=A[i][j]*B[j];
	  }
      }
    return ans;
  }

  template <typename T>
    TNT::Array1D<T> matmult(const TNT::Array1D<T>& A, const TNT::Array2D<T>& B)
  {//really using the transpose of A, produces a column vector when it 
    //should produce a row vector
    if(A.dim1()!=B.dim1())
      {
	throw invalid_argument("Mismatched matrix dimensions.");
      }
    TNT::Array1D<T> ans(B.dim2());
    
    for(int j=0; j<B.dim2(); j++)
      {
	ans[j]=0;
	for(int i=0; i<B.dim1(); i++)
	  {
	    ans[j]+=A[i]*B[i][j];
	  }
      }
    return ans;
  }


  template<typename T>
    void insert_1D(TNT::Array1D<T>& original, const TNT::Array1D<T>& inserted, const int origin)
    {
      if(origin+inserted.dim()>original.dim())
	{
	  throw out_of_range("Inserted array dimensions plus origin exceed \n original array dimensions.");
	}
      for(int j=origin; j<origin+inserted.dim(); j++)
	{
	  original[j]=inserted[j-origin];
	}
    }

  template <typename T>
    void insert_2D(TNT::Array2D<T>& original, const TNT::Array2D<T>& inserted, 
		   const int origin1, const int origin2)
    //origin gives the index of the upper left most element of the inserted array
    {
      if((origin1+inserted.dim1()>original.dim1())||
	 (origin2+inserted.dim2()>original.dim2()))
	{
	  throw out_of_range("Inserted array dimensions plus origin exceed \n original array dimensions.");
	}
      for(int i=origin1; i<origin1+inserted.dim1(); i++)
	{
	  for(int j=origin2; j<origin2+inserted.dim2(); j++)
	    {
	      original[i][j]=inserted[i-origin1][j-origin2];
	    }
	}
    }

  template <typename T>
    void insert_1D_into_2D(TNT::Array2D<T>& original, 
			   const TNT::Array1D<T>& inserted, 
			   const int index, bool iscolumn)
    //index gives the row or column into which the array should 
    //be inserted, iscolumn gives whether the 1D array is a column or row vector
    {
      if(iscolumn){
	if((index>=original.dim2())
	   ||(index<0)
	   ||(inserted.dim()!=original.dim1()))
	    {
	      throw out_of_range("Inserted array dimension not equal to external array dimension  or index out of range.");
	    }
	
	for(int i=0; i<inserted.dim(); i++)
	  {
	      original[i][index]=inserted[i];
	  }
	
      } else{//is row
	if((index>=original.dim1())
	   ||(index<0)
	   ||(inserted.dim()!=original.dim2()))
	  {
	    throw out_of_range("Inserted array dimension not equal to external array dimension or index out of range.");
	  }
	    for(int i=0; i<inserted.dim(); i++)
	      {
		original[index][i]=inserted[i];
	      }
      }    
    }

 
  template <typename T>
    void output1D(const TNT::Array1D<T>& matr)
    {
      for(int i=0;i<matr.dim();i++)
	{
	  cout << "( ";
	  cout << matr[i] <<" "; 
	  cout << ")\n";
	}
      cout << endl<< endl;
    }
  
  template <typename T>
    void output2D(const TNT::Array2D<T>& matr)
    {
      for(int i=0;i<matr.dim1();i++)
	{
	  cout << "( ";
      for(int j=0;j<matr.dim2();j++)
	{
	  cout << matr[i][j] <<" "; 
	}
      cout << ")\n";
	}
      cout << endl << endl;
    }

  
};


#endif 
