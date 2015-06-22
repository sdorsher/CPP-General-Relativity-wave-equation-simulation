#ifndef EVOLUTION_H
#define EVOLUTION_H

#include "TNT2.h"
#include "TwoDVectorGridFunction.h"
#include "GridFunction.h"
#include "Grid.h"
#include <fstream>
#include "ConfigParams.h"


//Low storage fourth order Runga Kutta routine.
//See pg 64 of Hesthaven and Warburton.

void rk4lowStorage(Grid thegrid, TwoDVectorGridFunction<double>& uh, 
                   TwoDVectorGridFunction<double>& RHSvgf, 
                   double t, double deltat);

#endif
