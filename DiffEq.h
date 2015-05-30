#ifndef DIFFEQ_H
#define DIFFEQ_H

#include "GridFunction.h"
#include <cmath>
#include "ConfigParams.h"


using namespace TNT;

GridFunction<Array2D<double>> setupAmatrix(GridFunction<double>& nodes);
GridFunction<Array2D<double>> setupBmatrix(GridFunction<double>& nodes);
#endif
