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
 public:
  double xp, xip, dxpdt, d2xpdt2, d2xdtdxi, Sminus, Splus;
  vector<bool> timeDepTrans;
  vector<double> dxdxib;
  vector<double> dxdxibL0;
  vector<double> dxdxibL1;
  vector<double> dxdxibR0;
  vector<double> dxdxibR1;
  
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
  

  void timedep_to_rstar(Orbit* orb);

  //HERE (OSCULATING ORBITS COORDINATE TRANSFORM)
  void coord_trans(double a, double b, Grid& thegrid, vector<double>&  x, vector<double> & dxdt, vector<double> & dxdxi, vector<double> & d2dxdt2, vector<double> & d2dxi2,vector<double> & d2xdtdxi, int elemnum);
  
    
};

#endif
