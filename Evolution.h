#include "tnt.h"
#include "TNT2.h"
#include "VectorGridFunction.h"
#include "GridFunction.h"

void rk4lowStorage(GridFunction& nodes, VectorGridFunction& uh, 
                   VectorGridFunction& RHSvgf, 
                   double t, double deltat);
void RHS(GridFunction& nodes, VectorGridFunction& uh, 
         VectorGridFunction& RHSvgf, double t);
