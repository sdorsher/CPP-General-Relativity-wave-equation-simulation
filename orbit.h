#ifndef ORBIT_H
#define ORBIT_H

#include "namespaces.h"
#include <cmath>
#include "ConfigParams.h"

using namespace orbit;

//set parameters in orbit namespace
void initialize_orbit();

//calculate phi as a function of time
double phi_of_t(double t);

//calculate chi as a function of time
double chi_of_t(double t);

#endif
