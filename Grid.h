#ifndef GRID_H
#define GRID_H

#include "tnt.h"
#include "TNT2.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include <fstream>

class Grid
{

 private:
  int NumElem;
  vector<double> elementBoundaries; //2N (should be N+1)?
  GridFunction nodeLocs; //NxNp
  //vector<Array2D> rescaledDmatrix;
  int order;
  vector<double> drdx; //jacobian
  //vector<int> elementOrders; //N
  void calcjacobian();

 public:
  Grid(string fileElemBoundaries, int elemorder, int numelems);
  //initializes elementBoundaries from file,
  //obtains node locations from reference element and puts them in array

  int numberElements();//return number of elements, calculated from input file
  GridFunction gridNodeLocations();  
  vector<double> gridBoundaries();
  Array2D<double> getRescaledDmatrix();
  vector<double> jacobian();

};


#endif