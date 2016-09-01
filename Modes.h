#ifndef MODES_H
#define MODES_H

#include <cmath>
#include <vector>
#include "namespaces.h"
#include "TwoDVectorGridFunction.h"
#include "gsl/gsl_sf_legendre.h"
#include "ConfigParams.h"
#include "Orbit.h"

//might require modifying to use Orbit.h instead of orbit.h


using namespace std;
using namespace orbit;

class Modes {
 public:
  vector<double> ll; //stores the l's to be used in computation
  vector<double> mm; //stores the corresponding m's to be used in computation
  vector<double> psil; // the scalar field, summed over m, by l mode
  vector<double> psitl; // the derivatie of the scalar field wrt t, summed over m, by l mode
  vector<double> psirl; // same thing wrt r
  vector<double> psiphil; // same thing wrt phi
  
  int ntotal; //total number of modes used in computation
  Modes(int lmax);
  void sum_m_modes(TwoDVectorGridFunction<complex<double>> uh,double time,int index1,int index2);
  // sum the modes to form psi, dpsi/dt, dpsi/dr, and dpsi/dphi at a certain radius
  
 private:
  int n_of_l(int l); //number of modes with l to be included in computation
  int nmodes_of_l(int lmax); //get ntotal
  void set_lm_mode_info(int lmax); //set l and m values, initialize other variables
};
#endif
