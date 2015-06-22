#ifndef GRID_H
#define GRID_H

#include "TNT2.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include <fstream>
#include "DiffEq.h"
#include "CharacteristicFlux.h"
#include "VectorGridFunction.h"

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
  GridFunction<Array2D<double>> Amatrices;
  VectorGridFunction<Array2D<double>> Bmatrices;
  int multipoleDim; 
  vector<double> ll; //mode list
  vector<CharacteristicFlux> AleftBoundaries;
  vector<CharacteristicFlux> ArightBoundaries;
  vector<Array1D<double>> duL; //A component of the flux, left boundary
  vector<Array1D<double>> duR; //A component of the flux, right boundary
  GridFunction<Array2D<double>> trimmedAmatrices;

  //Thoughts for the future:
  //vector<int> elementOrders; //N
  //map<int,ReferenceElement> allRefelems; 
  //to make multiple orders possible

 public:
  Grid(int elemorder, int numelements, int nummodes, double lowerlim, 
       double upperlim);
  ReferenceElement refelem; //member variable: the reference element

  int modesDim();
  int numberElements();//Returns number of elements, calculated from input file
  GridFunction<double> gridNodeLocations();  //Returns physical node location
  vector<double> gridBoundaries(); //Returns the boundaries of the elements
  double jacobian(int elemnum); //Returns the jacobian of a specific element

  //Returns du, for the characteristic flux.
  vector<TNT::Array2D<double>> characteristicflux(VectorGridFunction<double>& 
                                                  uh);
  //Returns the right hand side of the differential equation
  void RHS(VectorGridFunction<double>& uh, 
           VectorGridFunction<double>& RHSvgf, 
           double t, vector<Array2D<double>>& du );
};

#endif
