#ifndef GRID_H
#define GRID_H

#include "tnt.h"
#include "TNT2.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include <fstream>

class Grid
{

 public:
  int NumElem;
  vector<double> elementBoundaries; //2N (should be N+1)?
  GridFunction nodeLocations; //NxNp
  int order;
  
  //vector<int> elementOrders; //N

 public:
  Grid(string fileElemBoundaries, int elemorder, int numelems);
  //initializes elementBoundaries from file,
  //obtains node locations from reference element and puts them in array

  int getN();//return number of elements, calculated from input file
  GridFunction getNodeLocations();  
};


#endif
