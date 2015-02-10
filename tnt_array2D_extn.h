#ifndef TNT_ARRAY2D_EXTN_H
#define TNT_ARRAY2D_EXTN_H
#include <cmath>


#include "tnt_array2d.h"
//include "../tnt/tnt_array2d_utils.h"


namespace TNT
{


  template <class T>
    Array2D<T> sqrt(const Array2D<T> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++)
	{
	  for(int j=0; j<m; j++)
	    {
	      C[k][j]=std::sqrt(A[k][j]);
	    }
	}
      return C;
    }

  template <class T>
    Array2D<T> operator+(const Array2D<T> &A, const T B)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++)
	{
	  for(int j=0; j<m; j++)
	    {
	      C[k][j]=A[k][j]+B;
	    }
	}
      return C;
    }
  

  template <class T>
    Array2D<T> operator+(const T B, const Array2D<T> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++)
	{
	  for(int j=0; j<m; j++)
	    {
	      C[k][j]=A[k][j]+B;
	    }
	}
      return C;
    }


  template <class T>
    Array2D<T> operator*(const T B, const Array2D<T> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++)
	{
	  for(int j=0; j<m; j++)
	    {
	      C[k][j]=A[k][j]*B;
	    }
	}
      return C;
    }

  template <class T>
    Array2D<T> operator*(const Array2D<T> &A, const T B)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++)
	{
	  for(int j=0;j<m;j++)
	    {
	      C[k][j]=A[k][j]*B;
	    }
	}
      return C;
    }

  template <class T>
    Array2D<T> operator/(const Array2D<T> &A, const T B)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++)
	{
	  for(int j=0; j<m; j++)
	    {
	      C[k][j]=A[k][j]/B;
	    }
	}
      return C;
    }

  template <class T>
    Array2D<T> operator/(const T B, const Array2D<T> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++)
	{
	  for(int j=0; j<m; j++)
	    {
	      C[k][j]=B/A[k][j];
	    }
	}
      return C;
    }

  template <class T>
    Array2D<T> operator-(const T B, const Array2D<T> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++)
	{
	  for(int j=0; j<m; j++)
	    {
	      C[k][j]=B-A[k][j];
	    }
	}
      return C;
    }

  template <class T>
    Array2D<T> operator-(const Array2D<T> &A, const T B)
    {
      int n=A.dim1();
      int m=A.dim2();
	Array2D<T> C(n,m);
      for(int k=0; k<n; k++)
	{
	  for(int j=0; j<m; j++)
	    {
	      C[k][j]=A[k][j]-B;
	    }
	}
      return C;
    }

  

}



#endif
