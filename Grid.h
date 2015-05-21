#ifndef GRID_H
#define GRID_H

#include "TNT2.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include <fstream>
#include "DiffEq.h"
#include "CharacteristicFlux.h"

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
  vector<CharacteristicFlux> AleftBoundaries;
  vector<CharacteristicFlux> ArightBoundaries;
  //vector<int> elementOrders; //N
  //map<int,ReferenceElement> allRefelems;
  

 public:
  Grid(int elemorder, int numelements,double lowerlim, double upperlim);
  Grid(string fileElemBoundaries, int elemorder, int numelems);
  //initializes elementBoundaries from file,
  //obtains node locations from reference element and puts them in array

  int numberElements();//return number of elements, calculated from input file
  GridFunction<double> gridNodeLocations();  
  vector<double> gridBoundaries();
  Array2D<double> getRescaledDmatrix();
  vector<double> jacobian();
  ReferenceElement refelem;
};

#endif
