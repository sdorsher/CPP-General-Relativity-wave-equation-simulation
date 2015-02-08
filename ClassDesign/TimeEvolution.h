
#include "tnt_array1d.h"
#include "tnt_array2d."
#include "tnt_array1d_utils.h"
#include "tnt_array2d_utils.h"
#include "GridFunction.h"
#include "VectorGridFunction.h"
#include "DifferentialEquation.h"

class TimeEvolution
{
 private:
  Array2D coeffs;
  VectorGridFunction intermediateStep;
  DifferentialEquation de;
 public:
  void rk4lowStorage(GridFunction& nodes, VectorGridFunction& uh, 
		     VectorGridFunction& RHS);
  //calls de.RHS and de.boundaryConditions for each intermediate time step
};
    
