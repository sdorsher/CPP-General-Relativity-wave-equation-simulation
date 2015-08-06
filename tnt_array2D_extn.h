#ifndef TNT_ARRAY2D_EXTN_H
#define TNT_ARRAY2D_EXTN_H
#include <cmath>
#include <complex>

#include "tnt/tnt_array2d.h"
//include "../tnt/tnt_array2d_utils.h"

namespace TNT
{

  using namespace std;
  //Takes the transpose of the array
  template<class T>
    Array2D<T> transpose(const Array2D<T> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(m,n);
      for (int k=0; k<n; k++){
        for(int j=0; j<m; j++){
          C[j][k]=A[k][j];
        }
      }
      return C;
    }

  //Takes the square root of each element of the array
  template <class T>
    Array2D<T> sqrt(const Array2D<T> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++){
        for(int j=0; j<m; j++){
          C[k][j]=std::sqrt(A[k][j]);
        }
      }
      return C;
    }

  //Adds an array to a scalar
  template <class T>
    Array2D<T> operator+(const Array2D<T> &A, const T B)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++){
        for(int j=0; j<m; j++){
          C[k][j]=A[k][j]+B;
        }
      }
      return C;
    }
  
  //Adds a scalar to an array
  template <class T>
    Array2D<T> operator+(const T B, const Array2D<T> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<complex<T>> C(n,m);
      for(int k=0; k<n; k++){
        for(int j=0; j<m; j++){
          C[k][j]=A[k][j]+B;
        }
      }
      return C;
    }


  //Adds an array to a scalar
  template <class T>
    Array2D<complex<T>> operator+(const Array2D<complex<T>> &A, const T B)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<complex<T>> C(n,m);
      for(int k=0; k<n; k++){
        for(int j=0; j<m; j++){
          C[k][j]=A[k][j]+B;
        }
      }
      return C;
    }
  
  //Adds a scalar to an array
  template <class T>
    Array2D<complex<T>> operator+(const T B, const Array2D<complex<T>> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<complex<T>> C(n,m);
      for(int k=0; k<n; k++){
        for(int j=0; j<m; j++){
          C[k][j]=A[k][j]+B;
        }
      }
      return C;
    }
  //Multiplies a scalar by an array
  template <class T>
    Array2D<T> operator*(const T B, const Array2D<T> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++){
        for(int j=0; j<m; j++){
          C[k][j]=A[k][j]*B;
        }
      }
      return C;
    }


  //Multiplies an array by a scalar
  template <class T>
    Array2D<T> operator*(const Array2D<T> &A, const T B)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++){
	  for(int j=0;j<m;j++){
            C[k][j]=A[k][j]*B;
          }
      }
      return C;
    }
  //Multiplies a scalar by an array
  template <class T>
    Array2D<complex<T>> operator*(const T B, const Array2D<complex<T>> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<complex<T>> C(n,m);
      for(int k=0; k<n; k++){
        for(int j=0; j<m; j++){
          C[k][j]=A[k][j]*B;
        }
      }
      return C;
    }

  
  //Multiplies an array by a scalar
  template <class T>
    Array2D<complex<T>> operator*(const Array2D<complex<T>> &A, const T B)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++){
	  for(int j=0;j<m;j++){
            C[k][j]=A[k][j]*B;
          }
      }
      return C;
    }
  
  //Divides an array by a scalar
  template <class T>
    Array2D<T> operator/(const Array2D<T> &A, const T B)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++){
        for(int j=0; j<m; j++){
          C[k][j]=A[k][j]/B;
        }
      }
      return C;
    }

  //Divides a scalar by each element of an array
  template <class T>
    Array2D<T> operator/(const T B, const Array2D<T> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++){
        for(int j=0; j<m; j++){
          C[k][j]=B/A[k][j];
        }
      }
      return C;
    }

  //Subtracts an array from a scalar
  template <class T>
    Array2D<T> operator-(const T B, const Array2D<T> &A)
    {
      int n=A.dim1();
      int m=A.dim2();
      Array2D<T> C(n,m);
      for(int k=0; k<n; k++){
        for(int j=0; j<m; j++){
          C[k][j]=B-A[k][j];
        }
      }
      return C;
    }

  //Subtracts a scalar from an array
  template <class T>
    Array2D<T> operator-(const Array2D<T> &A, const T B)
    {
      int n=A.dim1();
      int m=A.dim2();
	Array2D<T> C(n,m);
      for(int k=0; k<n; k++){
        for(int j=0; j<m; j++){
          C[k][j]=A[k][j]-B;
        }
      }
      return C;
    }
}
#endif
