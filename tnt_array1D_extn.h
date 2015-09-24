#ifndef TNT_ARRAY1D_EXTN_H
#define TNT_ARRAY1D_EXTN_H
#include <cmath>
#include <complex>

#include "tnt/tnt_array1d.h"
//include "../tnt/tnt_array1d_utils.h"

namespace TNT
{

  using namespace std;
  //Takes the square root of each element
  template <class T>
    Array1D<T> sqrt(const Array1D<T> &A)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++){
        C[k]=std::sqrt(A[k]);
      }
      return C;
    }

  //Adds a sacalar to an array
  template <class T>
    Array1D<T> operator+(const Array1D<T> &A, const T B)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++){
        C[k]=A[k]+B;
      }
      return C;
    }
  
  //Adds a scalar to an array
  template <class T>
    Array1D<T> operator+(const T B, const Array1D<T> &A)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++){
        C[k]=A[k]+B;
      }
      return C;
    }
  //Adds a scalar to an array
  template <typename T>
    Array1D<complex<T>> operator+(const Array1D<complex<T>> &A, const T B)
    {
      int n=A.dim();
      Array1D<complex<T>> C(n);
      for(int k=0; k<n; k++){
        C[k]=A[k]+B;
      }
      return C;
    }
  
  //Adds a scalar to an array
  template <class T>
    Array1D<complex<T>> operator+(const T B, const Array1D<complex<T>> &A)
    {
      int n=A.dim();
      Array1D<complex<T>> C(n);
      for(int k=0; k<n; k++){
        C[k]=A[k]+B;
      }
      return C;
    }

  //Multiplies a scalar times an array
  template <class T>
    Array1D<T> operator*(const T B, const Array1D<T> &A)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++){
        C[k]=A[k]*B;
      }
      return C;
    }
  //Multiplies a scalar times a complex array
  template <class T>
    Array1D<complex<T>> operator*(const T B, const Array1D<complex<T>> &A)
    {
      int n=A.dim();
      Array1D<complex<T>> C(n);
      for(int k=0; k<n; k++){
        C[k]=A[k]*B;
      }
      return C;
    }

  //Multiplies a scalar times an array
  template <class T>
    Array1D<T> operator*(const Array1D<T> &A, const T B)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++){
        C[k]=A[k]*B;
      }
      return C;
    }
  
  //Divides an array by a scalar
  template <class T>
    Array1D<T> operator/(const Array1D<T> &A, const T B)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++){
        C[k]=A[k]/B;
      }
      return C;
    }

  //Divides a scalar by the elements of an array
  template <class T>
    Array1D<T> operator/(const T B, const Array1D<T> &A)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++){
        C[k]=B/A[k];
      }
      return C;
    }

  //Subtracts an array from a scalar
  template <class T>
    Array1D<T> operator-(const T B, const Array1D<T> &A)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++){
        C[k]=B-A[k];
      }
      return C;
    }

  //Subtracts a scalar from an array
  template <class T>
    Array1D<T> operator-(const Array1D<T> &A, const T B)
    {
      int n=A.dim();
      Array1D<T> C(n);
      for(int k=0; k<n; k++){
        C[k]=A[k]-B;
      }
      return C;
    }
}
#endif
