#ifndef MODES_H
#define MODES_H

#include <cmath>
#include <vector>
#include "namespaces.h"
#include "TwoDVectorGridFunction.h"
#include "gsl/gsl_sf_legendre.h"
#include "ConfigParams.h"

using namespace std;
using namespace orbit;

class Modes {
 public:
  vector<double> ll;
  vector<double> mm;
  vector<double> psil; // the scalar field, summed over m, by l mode
  vector<double> psitl; // the derivatie of the scalar field wrt t, summed over m, by l mode
  vector<double> psirl; // same thing wrt r
  vector<double> psiphil; // same thing wrt phi
  
  int ntotal;
  Modes(int lmax);
  void sum_m_modes(TwoDVectorGridFunction<complex<double>> uh,double time,int index1,int index2);

  
 private:
  int n_of_l(int l);
  int nmodes_of_l(int lmax);
  void set_lm_mode_info(int lmax);
};
#endif
