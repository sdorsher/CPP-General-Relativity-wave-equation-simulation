#ifndef TNT_ARRAY1D_EXTN_H
#define TNT_ARRAY1D_EXTN_H
#include <cmath>


#include "tnt/tnt_array1d.h"
//include "../tnt/tnt_array1d_utils.h"


namespace TNT
{

  template <class T>
    Array1D<T> sqrt(const Array1D<T> &A)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++)
	{
	  C[k]=std::sqrt(A[k]);
	}
      return C;
    }

  template <class T>
    Array1D<T> operator+(const Array1D<T> &A, const T B)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++)
	{
	  C[k]=A[k]+B;
	}
      return C;
    }
  

  template <class T>
    Array1D<T> operator+(const T B, const Array1D<T> &A)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++)
	{
	  C[k]=A[k]+B;
	}
      return C;
    }


  template <class T>
    Array1D<T> operator*(const T B, const Array1D<T> &A)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++)
	{
	  C[k]=A[k]*B;
	}
      return C;
    }

  template <class T>
    Array1D<T> operator*(const Array1D<T> &A, const T B)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++)
	{
	  C[k]=A[k]*B;
	}
      return C;
    }

  template <class T>
    Array1D<T> operator/(const Array1D<T> &A, const T B)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++)
	{
	  C[k]=A[k]/B;
	}
      return C;
    }

  template <class T>
    Array1D<T> operator/(const T B, const Array1D<T> &A)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++)
	{
	  C[k]=B/A[k];
	}
      return C;
    }

  template <class T>
    Array1D<T> operator-(const T B, const Array1D<T> &A)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++)
	{
	  C[k]=B-A[k];
	}
      return C;
    }

  template <class T>
    Array1D<T> operator-(const Array1D<T> &A, const T B)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++)
	{
	  C[k]=A[k]-B;
	}
      return C;
    }

  

}



#endif
