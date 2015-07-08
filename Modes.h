#ifndef MODES_H
#define MODES_H

#include <cmath>
#include <vector>

using namespace std;

class Modes {
 public:
  vector<double> ll;
  vector<double> mm;
  int ntotal;
  Modes(int lmax);
  
 private:
  int n_of_l(int l);
  int nmodes_of_l(int lmax);
  void set_lm_mode_info(int lmax);
};
#endif
