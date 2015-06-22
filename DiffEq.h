#ifndef DIFFEQ_H
#define DIFFEQ_H

#include "GridFunction.h"
#include <cmath>
#include "ConfigParams.h"

using namespace TNT;

//User supplied functions to set up the differential equation
// du/dt + A du/dx + Bu=0

GridFunction<Array2D<double>> setupABmatrices(GridFunction<double>& nodes);

#endif
