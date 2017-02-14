#ifndef GRID_H
#define GRID_H

#include "TNT2.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include <fstream>
#include "VectorGridFunction.h"
#include "TwoDVectorGridFunction.h"
#include <complex>
//#include "DiffEq.h"
#include "CharacteristicFlux.h"
#include "namespaces.h"
#include "OutputIndices.h"


using namespace std;
using namespace window;
using namespace layers;

class Grid
{
 private:
  int NumElem; //N. Number of elements in grid.
  int order;  //Np. Element order
  vector<double> elementBoundaries; //length N+1
  vector<double> drdx; //jacobian
  GridFunction<double> nodeLocs; //N by Np. Physical node locations.
  void calcjacobian();

 public:
  GridFunction<double> rschw;//Schwarzschild coordinate
  GridFunction<double> rstar; //tortoise coordinate

  GridFunction<double> window;
  GridFunction<double> dwindow;
  GridFunction<double> d2window;
  
  Grid(int elemorder, int numelements, int nummodes);
  ReferenceElement refelem; //member variable: the reference element

  //find the radii at which to output the data (returns grid and node indices)
  void find_extract_radii(double rfinite, double rSplus, OutputIndices& ijoutput, double dx);
  int numberElements();//Returns number of elements, calculated from input file
  GridFunction<double> gridNodeLocations();  //Returns physical node location
  vector<double> gridBoundaries(); //Returns the boundaries of the elements
  double jacobian(int elemnum); //Returns the jacobian of a specific element
  int nodeOrder(); //Returns the node order (same for all elements)
};

#endif
