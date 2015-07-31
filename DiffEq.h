#ifndef DIFFEQ_H
#define DIFFEQ_H

#include "source_interface.h"
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

using namespace TNT;
using namespace layers;
using namespace source_interface;

//User supplied functions to set up the differential equation
// du/dt + A du/dx + Bu=0

class DiffEq
{

 private:
  GridFunction<Array2D<double>> Amatrices;
  VectorGridFunction<Array2D<double>> Bmatrices;
  vector<CharacteristicFlux> AleftBoundaries;
  vector<CharacteristicFlux> ArightBoundaries;
  GridFunction<Array2D<double>> trimmedAmatrices;
  void setupABmatrices(Grid& thegrid, Modes& lmmodes);


 public:
  DiffEq(Grid& thegrid, Modes& lmmodes, int nmodetotal);
  Array2D<double> getA(int gridindex, int pointsindex);
  Array2D<double> getB(int modesindex, int gridindex, int pointsindex);
  CharacteristicFlux getAleft(int elemnum);
  CharacteristicFlux getAright(int elemnum);
  Array2D<double> getAtrimmed(int gridindex, int pointsindex);

  //Returns du, for the characteristic flux.
  vector<TNT::Array2D<complex<double>>> characteristicflux(int modenum, Grid& thegrid,
                                                  TwoDVectorGridFunction<complex<double>>& 
                                                  uh, bool output);
  //Returns the right hand side of the differential equation
  void RHS(int modenum, Grid& thegrid,
           TwoDVectorGridFunction<complex<double>>& uh, 
           TwoDVectorGridFunction<complex<double>>& RHSvgf, 
           double t, vector<Array2D<complex<double>>>& du , bool output);

  void modeRHS(Grid& thegrid,
               TwoDVectorGridFunction<complex<double>>& uh,
               TwoDVectorGridFunction<complex<double>>& RHStdgf, 
               double t, bool output); 
};

#endif
