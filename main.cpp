#include <iostream>
#include "TNT2.h"
#include "ReferenceElement.h"

using namespace TNT;

int main()
{
  int elemOrder = 5;
  ReferenceElement refelem(elemOrder);

  std::cout << "element order is " << refelem.getOrder() << std::endl;
  std::cout<< "reference node locations" <<std::endl;
  output1D(refelem.getr());
  std::cout << "-------------" << std::endl;
  std::cout << "-------------" << std::endl;
  std::cout << "derivative matrix" <<std::endl;
  output2D(refelem.getD());
  
}
