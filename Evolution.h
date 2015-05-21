#ifndef EVOLUTION_H
#define EVOLUTION_H

#include "TNT2.h"
#include "VectorGridFunction.h"
#include "GridFunction.h"
#include "Grid.h"
#include <fstream>
#include "ConfigParams.h"

void rk4lowStorage(Grid thegrid, VectorGridFunction<double>& uh, 
                   VectorGridFunction<double>& RHSvgf, 
                   double t, double deltat);

void RHS(Grid thegrid, VectorGridFunction<double>& uh, 
         VectorGridFunction<double>& RHSvgf, double t,bool output);

#endif
