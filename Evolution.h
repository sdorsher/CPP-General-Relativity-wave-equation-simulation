#ifndef EVOLUTION_H
#define EVOLUTION_H

#include "TNT2.h"
#include "VectorGridFunction.h"
#include "GridFunction.h"
#include "Grid.h"
#include <fstream>
#include "ConfigParams.h"


//low storage fourth order Runga Kutta routine
//see pg 64 of Hesthaven and Warburton

void rk4lowStorage(Grid thegrid, VectorGridFunction<double>& uh, 
                   VectorGridFunction<double>& RHSvgf, 
                   double t, double deltat);

#endif
