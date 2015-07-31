#ifndef ORBIT_H
#define ORBIT_H

#include "namespaces.h"
#include <cmath>
#include "ConfigParams.h"

using namespace orbit;

void initialize_orbit();
double phi_of_t(double t);
double chi_of_t(double t);

#endif
