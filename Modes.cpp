#include "Modes.h"

int n_of_l(int l)
{
  if (l % 2 == 0) {
    return (l + 2) / 2;
  } else {
    return (l + 1) / 2;
  }
}

int nmodes_of_l(int lmax)
{
  int nmodes = 0;
  
  for(int l=0; l<=lmax; l++) {
    nmodes += n_of_l(l);
  }
  return nmodes;
}

  
