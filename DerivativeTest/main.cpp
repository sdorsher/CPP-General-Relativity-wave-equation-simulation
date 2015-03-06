#include "../TNT2.h"
#include "../ReferenceElement.h"
#include <cmath>
#include <iostream>
#include <fstream>

//test the reference element derivative matrix on a sinusoidal function,
//plot convergence of the L_infinity error with order of the element
// L_infinity error is the max error

using namespace TNT;

int main()
{
  double pi = 3.14159265359;
  double period = 2.0; //-1 to 1 (from the ends of the reference element)
  double frequency = 1.0/period;
  double phase = 0.0;
  double amplitude =1.0;
  double dt = 0.01;
  int maxorder = 99;
  int minorder = 6;
  vector<ReferenceElement> elems;

  ofstream fs;
  fs.open("sinusoidDerivative.txt");
  
  for(int order=minorder; order<=maxorder; order++)
    {
      ReferenceElement refelem(order);
      elems.push_back(refelem);
      Array1D<double> r = elems[order-minorder].getr();
      Array1D<double> sinusoid(order+1);
      Array1D<double> cosinisoid(order+1);
      for(int i=0; i<=order; i++)
	{
	  sinusoid[i]= amplitude*sin(2.0*pi*frequency*r[i]+phase);
	  cosinisoid[i]= amplitude*2.0*pi*frequency
	    *cos(2.0*pi*frequency*r[i]+phase);
	}
	Array2D<double> Dmatrix(order,order); 
	Dmatrix = elems[order-minorder].getD();
      Array1D<double> Dsinusoid= matmult(Dmatrix,sinusoid);
      double maxerror=0.0;
      for(int i = 0; i<=order; i++)
	{
	  double errori = abs(cosinisoid[i] - Dsinusoid[i]);
	  maxerror = (errori >= maxerror) ? errori : maxerror;
	}
      fs << order << "\t" << maxerror << "\n";
    }
  fs.close();
 }
