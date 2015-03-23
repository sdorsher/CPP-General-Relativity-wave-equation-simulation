#include <iostream>
#include "TNT2.h"
#include "ReferenceElement.h"
#include "GridFunction.h"

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
  
  Array1D<double> threes(3,3.0);
  Array1D<double> twos(2,2.0);
  GridFunction gf(8,2,true);
  output1D(gf.get(3));
  gf.set(3,twos);
  output1D(gf.get(3));
  gf.set(3,0,3.0);
  output1D(gf.get(3));  
  //  gf.set(3,threes);
  //  gf.set(8,twos);
  //gf.set(8,0,3.0);
  //gf.set(3,3,3.0);

  gf.append(twos);
  output1D(gf.get(7));
  output1D(gf.get(8));

  cout << gf.get(3,0) << endl;
  //  gf.get(9,0);
  //gf.get(3,3);

  cout << gf.gridDim() << "\t" << gf.pointsDim() << endl;

  output1D(gf.get(8));

  cout << "-------------------" << endl;
  GridFunction gf2(0,2,false);
  gf2.append(twos);
  output1D(gf2.get(0));
  cout << "------------------"<<endl;
  GridFunction gfnew = gf+gf;
  output1D(gfnew.get(3));

  gfnew.save("testgf.txt");
  
}
