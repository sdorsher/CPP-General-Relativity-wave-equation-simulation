#ifndef EVOLUTION_H
#define EVOLUTION_H

#include "TNT2.h"
#include "TwoDVectorGridFunction.h"
#include "GridFunction.h"
#include "DiffEq.h"
#include <fstream>
#include "ConfigParams.h"
#include <complex>

//Low storage fourth order Runga Kutta routine.
//See pg 64 of Hesthaven and Warburton.

using namespace std;

//time evolution: Fourth order runga kutta 
void rk4lowStorage(Grid thegrid, DiffEq theequation, 
                   TwoDVectorGridFunction<complex<double>>& uh, 
                   TwoDVectorGridFunction<complex<double>>& RHSvgf, 
                   double t, double deltat);

#endif
