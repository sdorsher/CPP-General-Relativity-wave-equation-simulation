#ifndef HYPERBOLOIDALCOORDS_H
#define HYPERBOLOIDALCOORDS_H

#include "globals.h"
#include <cmath>
#include <iostream>

using namespace std;

double Lambert(double z);
double rschw(double z, double mass);
double invert_tortoise(double rstar, double mass);
void transition(double rho, double R, double S, double& fT, double& fTp,
                double& fTpp);

#endif
