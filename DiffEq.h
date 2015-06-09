#ifndef DIFFEQ_H
#define DIFFEQ_H

#include "GridFunction.h"
#include <cmath>
#include "ConfigParams.h"

using namespace TNT;

//user supplied functions to set up the differential equation
// du/dt + A du/dx + Bu=0

GridFunction<Array2D<double>> setupAmatrix(GridFunction<double>& nodes);
GridFunction<Array2D<double>> setupBmatrix(GridFunction<double>& nodes);

#endif
