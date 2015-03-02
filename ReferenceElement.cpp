#include "ReferenceElement.h"
#include <iostream>

using namespace TNT;

ReferenceElement::ReferenceElement(int N):refNodeLocations(N+1),vandermondeMatrix(N+1,N+1),dVdr(N+1,N+1),
					  derivativeMatrix(N+1,N+1)
{
  order=N;
 refNodeLocations=jacobiGL(0.0,0.0,N);
 vandermonde1D();
  gradVandermonde1D();
  Dmatrix1D();
  
}
					   //{

//  refNodeLocations=jglQuadraturePoints();
//  vandermondeMatrix=computeVandermonde();
//  dVdr=computedVdr();
//  derivativeMatrix=computeDr();
					   //}

void ReferenceElement::jacobiGQ(Array1D<double>& x, double alpha, 
				double beta, int n, Array1D<double>& w)
//checks out for n=5,7,9
  
{

  if((w.dim()!=n+1)||x.dim()!=n+1)
    {
      throw std::invalid_argument("JacobiGQ: Array1D argument dimensions do not agree with argument order n.\n");
    }

  Array1D<double> h1(n+1);
  Array1D<double> jdiag1(n+1);
  Array1D<double> idouble(n+1);
  Array1D<double> jdiag2(n);
  Array1D<double> d1(n);
  Array1D<double> d2(n);
  Array1D<double> d3(n);
  Array2D<double> vect(n+1,n+1); //contains eigenvectors
  double lngammaab, lngammaa, lngammab;

  if (n==0)
    {
      x[0]=(alpha-beta)/(alpha+beta+2.0);
      w[0]=2.0;
      return;
    }

  for(int i=0; i<=n; i++)
    {
      idouble[i] = i*1.0;
    }

  h1=2.0*idouble+alpha+beta;
  for(int i=0; i<=n; i++)
    {
      if(h1[i]>0.0)
	{
	  jdiag1[i]=-(pow(alpha,2.0)-pow(beta,2.0))/(h1[i]*(h1[i]+2.0));
	}				   
      else 
	{
	  jdiag1[i]=0.0;
	}
    }

  for(int i=0; i<n; i++)
    {
      d1[i]=2.0/(h1[i]+2.0);
      d2[i]=idouble[i+1]*(idouble[i+1]+alpha+beta)*(idouble[i+1]+alpha)*(idouble[i+1]+beta);
      d3[i]=1.0/((h1[i]+1.0)*(h1[i]+3.0));
    }
  
  jdiag2=d1*sqrt(d2*d3);

  //tridiagonal matrix has jdiag1 along the diagonal and jdiag2 to either side, zeros elsewhere

  Array2D<double> tridiag(n+1,n+1,0.0);

  for(int i=0; i<=n; i++)
    {
      for(int j=0; j<=n; j++)
	{
	  if(i==j)
	    {
	      tridiag[i][j]=jdiag1[i];
	    }
	  else if(abs(i-j)==1)
	    {
	      tridiag[i][j]=jdiag2[std::min(i,j)];
	    }
	}
    }
  //vect eigenvector, jdiag1 eigenvalues

  JAMA::Eigenvalue<double> eigen(tridiag);
  eigen.getRealEigenvalues(x);

  eigen.getV(vect);
  lngammaab=lgamma(alpha+beta+1.0);
  lngammaa=lgamma(alpha+1.0);
  lngammab=lgamma(beta+1.0);
  for(int i=0; i<=n; i++)
    {
      w[i]=pow(vect[0][i],2.0)*pow(2.0,(alpha+beta+1.0))
	/(alpha+beta+1.0)*exp(lngammaa+lngammab-lngammaab);
    }

  
}

Array1D<double> ReferenceElement::jacobiGL(double alpha, double beta, double n)
{
  //sets Jacobi-Gauss-Lebato quadrature points

  Array1D<double> quadpts(n+1,0.0);
  Array1D<double> w(n-1); //weights
  Array1D<double> x(n-1); //JacobiGQ input and output

  if(n<=0)
    throw std::invalid_argument("JacobiGL called with n<=0. Aborting");
  
  if(n==1)
    {
      quadpts[0]=-1.0;
      quadpts[1]=1.0;
      return quadpts;
    }else
    {
      jacobiGQ(x, alpha+1.0, beta+1.0, n-2, w);
      insert_1D(quadpts,x,1);
      quadpts[0]=-1.0;
      quadpts[n]=1.0;
      return quadpts;
    }
}


Array1D<double> ReferenceElement::jacobiP(const Array1D<double>& x, double alpha,  
			       double beta, int n)
//Evaluate Jacobi Polynomial of type (alpha,beta) at points x for order N
//may be able to store memory by not storing pl
// tested for alpha=beta=0, sizex=4 (n=0 through 3)
// tested for alpha=beta=1 sizex=4 (n=0 through 3)


{

  Array1D<double> polynom(x.dim());
  vector<Array1D<double>> pl;
  double lngammaab=lgamma(alpha+beta+1.0);
  double lngammaa =lgamma(alpha+1.0);
  double lngammab=lgamma(beta+1.0);

  pl.resize(n+1);

  double invsqgamma0=pow(2.0,alpha+beta+1.0)/(alpha+beta+1.0)
    *exp(lngammaa+lngammab-lngammaab);

  double gamma0=1.0/sqrt(invsqgamma0);

  Array1D<double> gamma0arr(x.dim(),gamma0);

  if (n==0)
    {
      polynom= gamma0arr;
      return polynom;
    }
  pl[0]=gamma0arr;
  

  double gamma1 = 0.5*sqrt((alpha+beta+3.0)/((alpha+1.0)*(beta+1.0)))*gamma0;
  double fac1 = (alpha+beta+2.0);
  double fac2 = (alpha-beta);
  //  std::cout << "gamma1 " << gamma1 << " fac1 " << fac1 << " fac2 " << fac2 <<std::endl;
  //checks out up to here
  
  Array1D<double> gamma1arr= gamma1*(fac1*x+fac2);


  if (n==1)
    {
      polynom= gamma1arr;
      return polynom;
    }	    
  
  pl[1]=gamma1arr;

  double aold = 2.0/(2.0+alpha+beta)*sqrt((1.0+alpha)*(1.0+beta)
					 /(3.0+alpha+beta));
  for(int i=0;i<=n-2;i++)
    {
      double idouble= double(i)+1.0;
      double idoublep1=idouble+1.0;
      double h1 = 2.0*idouble+alpha+beta;
      double anew = 2.0/(h1+2.0)*sqrt(idoublep1*(idoublep1+alpha+beta)*
				      (idoublep1+alpha)*(idoublep1+beta)
				      /(h1+1.0)/(h1+3.0));
      double bnew = -(alpha*alpha-beta*beta)/(h1*(h1+2.0));
      pl[i+2]=1.0/anew*(-aold* pl[i]+(x-bnew)*pl[i+1]);

      aold=anew;
 
    }//end for
  
  polynom= pl[n];
  return polynom;

}

void ReferenceElement::vandermonde1D()
{
  for(int j = 0; j<=order; j++)
    {
      insert_1D_into_2D(vandermondeMatrix,jacobiP(refNodeLocations,0.0,0.0,j),j,true);
    }


}


Array1D<double> ReferenceElement::gradJacobiP(double alpha, double beta,int N)
//evaluated the derivative of the Jacobi polynomials of type alpha, beta >-1 at the nodes for order N 

{
  if ((alpha<-1)||(beta<-1)) throw invalid_argument("alpha or beta <-1 in gradJacobiP.");
  if (N<0) throw std::invalid_argument("N<0 in gradJacobiP");
  
  Array1D<double> gradpoly(refNodeLocations.dim(),0.0);
  if (N!=0)
    {
      gradpoly = sqrt(N*(N+alpha+beta+1))*jacobiP(refNodeLocations,alpha+1.0,beta+1.0,N-1);
    }
  return gradpoly;
}

void ReferenceElement::gradVandermonde1D()
{
 
  for(int i=0; i<=order;i++)
    {
      insert_1D_into_2D(dVdr,gradJacobiP(0.0,0.0,i),i,true);
    }
  return;

}


void ReferenceElement::Dmatrix1D()
{

  Array2D<double> VT = transpose(vandermondeMatrix);
JAMA::LU<double> solver(VT);
  Array2D<double> dVdrT = transpose(dVdr);
  Array2D<double> DT = solver.solve(dVdrT);
  derivativeMatrix = transpose(DT);

//  JAMA::LU<double> solver(transpose(vandermondeMatrix);
//  derivativeMatrix=transpose(solver.solve(transpose(dVdr)));

}
