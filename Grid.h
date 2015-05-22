#ifndef GRID_H
#define GRID_H

#include "TNT2.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include <fstream>
#include "DiffEq.h"
#include "CharacteristicFlux.h"
#include "VectorGridFunction.h"

class Grid
{
 private:
  int NumElem;
  int order;  
  vector<double> elementBoundaries; //(N+1)
  vector<double> drdx; //jacobian
  GridFunction<double> nodeLocs; //NxNp
  void calcjacobian();
  GridFunction<Array2D<double>> Amatrices;
  GridFunction<Array2D<double>> Bmatrices;
  vector<CharacteristicFlux> AleftBoundaries;
  vector<CharacteristicFlux> ArightBoundaries;
  //vector<int> elementOrders; //N
  //map<int,ReferenceElement> allRefelems;
  vector<Array1D<double>> duL;
  vector<Array1D<double>> duR;



 public:
  Grid(int elemorder, int numelements,double lowerlim, double upperlim);
  Grid(string fileElemBoundaries, int elemorder, int numelems);
  //initializes elementBoundaries from file,
  //obtains node locations from reference element and puts them in array

  int numberElements();//return number of elements, calculated from input file
  GridFunction<double> gridNodeLocations();  
  vector<double> gridBoundaries();
  Array2D<double> getRescaledDmatrix();
  double jacobian(int elemnum);
  ReferenceElement refelem;
  vector<TNT::Array2D<double>> characteristicflux(VectorGridFunction<double>& uh);
  void RHS(VectorGridFunction<double>& uh, 
           VectorGridFunction<double>& RHSvgf, double t, vector<Array2D<double>>& du );
  

};

#endif
