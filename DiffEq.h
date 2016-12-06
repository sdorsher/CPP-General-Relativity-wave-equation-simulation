#ifndef DIFFEQ_H
#define DIFFEQ_H

#include "GridFunction.h"
#include <cmath>
#include "ConfigParams.h"
#include "CharacteristicFlux.h"
#include "Grid.h"
#include "VectorGridFunction.h"
#include "Modes.h"
#include <cfloat>
#include "namespaces.h"
#include <complex>
#include "source_interface.h"
#include <iomanip>
#include <omp.h>
#include "vecMatrixTools.h"
#include "TwoDVectorGridFunction.h"
#include "Coordinates.h"
#include "Orbit.h"
#include "EllipticalOrbit.h"
#include "CircularOrbit.h"
#include "WorldTube.h"

using namespace TNT;
using namespace layers;
using namespace source_interface;
using namespace std;
using namespace orbit;

//User supplied functions to set up the differential equation
// du/dt + A du/dx + Bu=0
// A is trimmed by cutting out the rows that are all zero. It should be square. 


class DiffEq
{

 private:
  //Grid functions and VectorGridFunctions range over all elements and all nodes
  VectorGridFunction<double> Amatrices; //used to be GridFunction array2D
  TwoDVectorGridFunction<double> Bmatrices; //used to be Vector Grid Function array2D
  VectorGridFunction<double> trimmedAmatrices; //used to be GridFunction array2D
  vector<complex<double>> uintL; //internal u at left boundary
  vector<complex<double>> uintR; //internal u at right boundary
  vector<complex<double>> uextL; //external u at left boundary
  vector<complex<double>> uextR; //external u at right boundary
  //vectors range over the elements, containing left and right boundaries
  vector<CharacteristicFlux> AleftBoundaries;
  vector<CharacteristicFlux> ArightBoundaries;

  //the function to setup the A and B matricies, written specificially for
  //the differential equation in question (with options based on parameters
  //for different differential equations)
  void setupABmatrices(Grid& thegrid, Modes& lmmodes, Coordinates& coordobj);

 public:
  VectorGridFunction<complex<double>> source; //the effective source
  
 public:
  DiffEq(Grid& thegrid, Modes& lmmodes, int nmodetotal, Coordinates & coords);
  vector<double> getA(int gridindex, int pointsindex);
  vector<double> getB(int modesindex, int gridindex, int pointsindex);
  CharacteristicFlux getAleft(int elemnum);
  CharacteristicFlux getAright(int elemnum);
  vector<double> getAtrimmed(int gridindex, int pointsindex);

  void set_coefficients(Grid &thegrid, EllipticalOrbit* orb, Coordinates &coords, double& maxspeed, int elemnum, Modes& lmmodes, double xp, double xip, double dxpdt, double d2xpdt2, vector<double>& dxdt, vector<double>& dxdxi);
  
  //Returns du, for the characteristic flux.
  vector<vector<complex<double>>> characteristicflux(double t,
							   int modenum,
							   Grid& thegrid,
				   TwoDVectorGridFunction<complex<double>>& 
							   uh, bool output);
  //Returns the right hand side of the differential equation including the flux
  void RHS(int modenum, Grid& thegrid,
           TwoDVectorGridFunction<complex<double>>& uh, 
           TwoDVectorGridFunction<complex<double>>& RHStdvgf, 
           double t, bool output, Coordinates& coords, WorldTube* wt);

  //loops over the modes to get the effective source,
  //get the characteristic flux, then calculate the RHS of the
  //differential equation
  void modeRHS(Grid& thegrid,
               TwoDVectorGridFunction<complex<double>>& uh,
               TwoDVectorGridFunction<complex<double>>& RHStdgf, 
               double t, bool output, Orbit* orb, WorldTube* wt, Coordinates& coords, double & max_speed, Modes& lmmodes); 

};

#endif
