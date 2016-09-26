#ifndef EVOLUTION_H
#define EVOLUTION_H

#include "TwoDVectorGridFunction.h"
#include "GridFunction.h"
#include "DiffEq.h"
#include <fstream>
#include "ConfigParams.h"
#include <complex>
#include "source_interface.h"
#include "WorldTube.h"
#include "Orbit.h"
#include "EllipticalOrbit.h"
#include "CircularOrbit.h"
#include "Coordinates.h"

//Low storage fourth order Runga Kutta routine.
//See pg 64 of Hesthaven and Warburton.

using namespace std;

//time evolution: Fourth order runga kutta 
void rk4lowStorage(Grid thegrid, DiffEq theequation, 
                   TwoDVectorGridFunction<complex<double>>& uh, 
                   TwoDVectorGridFunction<complex<double>>& RHSvgf, 
                   double t, double deltat, WorldTube* wt, double & max_speed, Orbit * orb, Coordinates & coords, Modes& lmmodes);

#endif
