#ifndef HYPERBOLOIDALCOORDS_H
#define HYPERBOLOIDALCOORDS_H

#include "globals.h"
#include <cmath>
#include <iostream>

using namespace std;

double rstar_of_r(double r,double mass); //calculate tortoise coordinate from schwarszchild coordinates
double Lambert(double z);// Function to calculate Lambert's W-function
double rschw(double z, double mass);//Inverts tortoise radius to get rschw-2M using Newton root finding method
double invert_tortoise(double rstar, double mass);//Inverts tortiose radius to get rschw-2M, using method dependent upon regime
//the transition between the hyperboloidal and tortoise regions
void transition(double rho, double R, double S, double& fT, double& fTp,
                double& fTpp);


#endif
