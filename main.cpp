#include <iostream>
#include "TNT2.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include "VectorGridFunction.h"
#include "Grid.h"

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
 
  cout << "------------------"<<endl;
  cout << "------------------"<<endl;
  cout << "begin vector grid function" <<endl;

  VectorGridFunction vgf(3,8,2,true);
  output1D(vgf.data.at(1).get(3));
  
  cout << vgf.vectorDim() << "\t" << vgf.gridDim() << "\t" << vgf.pointsDim() <<endl;

  vgf.set(0,gf);
  vgf.set(1,0,twos);
  vgf.set(2,0,0,3.0);
  GridFunction gftest=vgf.get(2);
  output1D(gftest.get(0));
  output1D(vgf.get(0,3));
  cout << vgf.get(1,0,0) <<endl;
  //successfull

  VectorGridFunction vgfnew = vgf+vgf;
  output1D(vgfnew.get(0,3));

  vgfnew.save("testvgf.txt");

  vector<double> outputvec = vgfnew.getVector(0,0);
  for(int i=0; i<8; i++)
    {
      cout << outputvec[i] <<endl;
    }

  void (*pf)(GridFunction&, VectorGridFunction&, VectorGridFunction&, int, int);
  pf=&testfunc;
  Array1D<double> ones(2,1.0);
  GridFunction gfones(0,2,false);
  for (int i=0; i<5; i++)
    {
      gfones.append(ones);
    }
  VectorGridFunction uh(3,5,2,true);
  VectorGridFunction RHS(3,5,2,true);
  loop(gfones,uh,RHS,pf);
  RHS.save("RHStest.txt");

  cout << "-------------"<<endl;
  cout << "begin grid test"<<endl;
  
  Grid gr("elemBoundaries.txt",20,13);
  vector<double> elemBounds = gr.gridBoundaries();
  for(int elem=0; elem<=gr.numberElements(); elem++)
    {
      cout << elemBounds[elem] << endl;
    }

  gr.gridNodeLocations().save("physicalNodes.txt");
  //Grid works

  
}
