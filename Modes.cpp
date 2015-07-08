#include "Modes.h"

Modes::Modes(int lmax) 
{
  set_lm_mode_info(lmax);
}

int Modes::n_of_l(int l)
{
  if (l % 2 == 0) {
    return (l + 2) / 2;
  } else {
    return (l + 1) / 2;
  }
}

int Modes::nmodes_of_l(int lmax)
{
  int nmodes = 0;
  
  for(int l=0; l<=lmax; l++) {
    nmodes += n_of_l(l);
  }
  return nmodes;
}

void Modes::set_lm_mode_info(int lmax) {
  ntotal = nmodes_of_l(lmax);
  ll.resize(ntotal);
  mm.resize(ntotal);
  int n = 0;
  for(int l=0; l<=lmax; l++){
    for(int m = l%2; m<=l; m+=2){
      ll[n]=l;
      mm[n]=m;
      n++;
    }
  }
}
