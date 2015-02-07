#include <iostream>
#include "../tnt2/TNT2.h"
#include "ReferenceElement.h"

using namespace TNT;

int main()
{
  int elemOrder = 5;
  ReferenceElement refelem(elemOrder);
  
  Array1D<double> x(elemOrder+1);
  double elemLeft=-1.0;
  double elemRight=1.0;
  double step=(elemRight-elemLeft)/double(elemOrder);
 
  double alpha =0.0;
  double beta =0.0;
  int n=5;

  Array1D<double> w(elemOrder+1);
  
  //Test JacobiGL
  
  std::cout << "test JacobiGL" << std::endl;
  x=refelem.jacobiGL(alpha,beta,n);
  output1D(x);
  
  std::cout << "-------------" << std::endl;
  std::cout << "-------------" << std::endl;


  //Test JacobiGQ-- n always equals elemOrder-1
  std::cout << "test JacobiGQ" << std::endl;
  refelem.jacobiGQ(x, alpha, beta, n, w);

  std::cout << "x:" << std::endl;
  output1D(x);
  std::cout << "-------" << std::endl;
  std::cout << "w:" << std::endl;
  output1D(w);
 

  
  std::cout<< "------------" <<std::endl;
  std::cout<< "------------" <<std::endl;
  std::cout << "test JacobiP" << std::endl;

  n=3;

    //test JacobiP 
    for(int i=0;i<=elemOrder+1;i++)
    {
      x[i]=elemLeft+step*double(i);
    }
  
  Array1D<double> polynom(elemOrder+1);


  refelem.jacobiP(x,alpha,beta,n,polynom);
  output1D(polynom);
  
}
