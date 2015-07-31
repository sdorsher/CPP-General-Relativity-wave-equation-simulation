#ifndef GRID_H
#define GRID_H

#include "TNT2.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include <fstream>
#include "VectorGridFunction.h"
#include "TwoDVectorGridFunction.h"
#include <complex>

// du/dt + A du/dx + Bu = 0
// See DiffDeq.cpp for A,B definition
// A is trimmed by cutting out the rows that are all zero.

class Grid
{
 private:
  int NumElem; //N. Number of elements in grid.
  int order;  //Np. Element order
  vector<double> elementBoundaries; //(N+1)
  vector<double> drdx; //jacobian
  GridFunction<double> nodeLocs; //N by Np. Physical node locations.
  void calcjacobian();

  //Thoughts for the future:
  //vector<int> elementOrders; //N
  //map<int,ReferenceElement> allRefelems; 
  //to make multiple orders possible

 public:
  GridFunction<double> rschw;
  GridFunction<double> rstar;
  VectorGridFunction<complex<double>> source;
  GridFunction<double> window;
  GridFunction<double> dwindow;
  GridFunction<double> d2window;
  Grid(int elemorder, int numelements, int nummodes, double lowerlim, 
       double upperlim);
  ReferenceElement refelem; //member variable: the reference element

  void find_extract_radii(double rfinite, double rSplus, int& ifinite, 
                              int& iSplus, int& jfinite, int& jSplus);
  int numberElements();//Returns number of elements, calculated from input file
  GridFunction<double> gridNodeLocations();  //Returns physical node location
  vector<double> gridBoundaries(); //Returns the boundaries of the elements
  double jacobian(int elemnum); //Returns the jacobian of a specific element
  int nodeOrder(); //Returns the node order (same for all elements)

};

#endif
