#ifndef DIFFEQ_H
#define DIFFEQ_H

#include "GridFunction.h"
#include <cmath>
#include "ConfigParams.h"
#include "CharacteristicFlux.h"
#include "Grid.h"
#include "VectorGridFunction.h"
#include "Modes.h"
#include "HyperboloidalCoords.h"
#include <cfloat>
#include "namespaces.h"
#include <complex>
#include "source_interface.h"
#include <iomanip>
#include <omp.h>

using namespace TNT;
using namespace layers;
using namespace source_interface;
using namespace std;

//User supplied functions to set up the differential equation
// du/dt + A du/dx + Bu=0
// A is trimmed by cutting out the rows that are all zero. It should be square. 


class DiffEq
{

 private:
  //Grid functions and VectorGridFunctions range over all elements and all nodes
  GridFunction<Array2D<double>> Amatrices;
  VectorGridFunction<Array2D<double>> Bmatrices;
  GridFunction<Array2D<double>> trimmedAmatrices;

  //vectors range over the elements, containing left and right boundaries
  vector<CharacteristicFlux> AleftBoundaries;
  vector<CharacteristicFlux> ArightBoundaries;

  //the function to setup the A and B matricies, written specificially for
  //the differential equation in question (with options based on parameters
  //for different differential equations)
  void setupABmatrices(Grid& thegrid, Modes& lmmodes);

 public:
  VectorGridFunction<complex<double>> source; //the effective source
  GridFunction<double> window; //the window function that factors into the effective source
  GridFunction<double> dwindow; //derivative of the window function
  GridFunction<double> d2window; //second derivative of the window function

 public:
  DiffEq(Grid& thegrid, Modes& lmmodes, int nmodetotal);
  Array2D<double> getA(int gridindex, int pointsindex);
  Array2D<double> getB(int modesindex, int gridindex, int pointsindex);
  CharacteristicFlux getAleft(int elemnum);
  CharacteristicFlux getAright(int elemnum);
  Array2D<double> getAtrimmed(int gridindex, int pointsindex);

  //Returns du, for the characteristic flux.
  vector<TNT::Array2D<complex<double>>> characteristicflux(double t,
							   int modenum,
							   Grid& thegrid,
				   TwoDVectorGridFunction<complex<double>>& 
							   uh, bool output);
  //Returns the right hand side of the differential equation including the flux
  void RHS(int modenum, Grid& thegrid,
           TwoDVectorGridFunction<complex<double>>& uh, 
           TwoDVectorGridFunction<complex<double>>& RHSvgf, 
           double t, vector<Array2D<complex<double>>>& du , bool output);

  //loops over the modes to get the effective source,
  //get the characteristic flux, then calculate the RHS of the
  //differential equation
  void modeRHS(Grid& thegrid,
               TwoDVectorGridFunction<complex<double>>& uh,
               TwoDVectorGridFunction<complex<double>>& RHStdgf, 
               double t, bool output); 
};

#endif
