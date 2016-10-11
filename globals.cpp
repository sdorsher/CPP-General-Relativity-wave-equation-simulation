#include "globals.h"

double PI=acos(-1.0);

int SOURCE_VECNUM = 2; //vector index that the effective source effects on the
//right hand side of the differential equation. Corresponds to time derivative.
//params is also a global, located in ConfigParams.h

//Indices of the grid (i) and of the nodes (j) where data will be output
// at a finite radius and at Splus. Based on params.grid.outputradius

