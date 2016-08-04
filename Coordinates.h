#ifndef COORDINATES_H
#define COORDINATES_H

#include "globals.h"
#include <cmath>
#include <iostream>

using namespace std;

class Coordinates{

  double R1, R2, W1, W2;
  
  Coordinates();
  
  double rstar_of_r(double r,double mass); //calculate tortoise coordinate from
  //Schwarszchild coordinates
  
  double Lambert(double z);// Function to calculate Lambert's W-function
  
  double rschw(double z, double mass);//Inverts tortoise radius to get rschw-2M using
  //Newton root finding method

  double invert_tortoise(double rstar, double mass);//Inverts tortiose radius to get rschw-2M,
  //using method dependent upon regime
  
  void transition(double rho, double R, double S, double& fT, double& fTp,
		  double& fTpp);
  //the transition between the hyperboloidal and tortoise regions
  
  
}

#endif
