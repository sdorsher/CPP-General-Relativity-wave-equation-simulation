#ifndef COORDINATES_H
#define COORDINATES_H

#include "globals.h"
#include <cmath>
#include <iostream>
#include "ConfigParams.h"
#include <vector>
#include "Orbit.h"
#include "CircularOrbit.h"
#include "EllipticalOrbit.h"
#include "Grid.h"


using namespace std;

class Coordinates{

  double R1, R2, W1, W2, a, b, xp, xip, dxpdt, d2xpdt2, d2xdtdxi;
  vector<double> dxdxib;
  vector<double> dxidbL;
  vector<double> dxidbR;
  
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
  

  void timedep_to_rstar(Orbit& orb);

  //HERE (OSCULATING ORBITS COORDINATE TRANSFORM)
  void coord_trans(Coordinates &coords, Grid& thegrid, vector<double> & dxdxi, vector<double> & d2dxdt2, vector<double> & d2dxdxi2,vector<double> & d2xdtdxi, int elemnum);
  
    
};

#endif
