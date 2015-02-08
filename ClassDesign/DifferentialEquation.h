#include "tnt_array1d.h"
#include "tnt_array2d."
#include "tnt_array1d_utils.h"
#include "tnt_array2d_utils.h"

class DifferentialEquation
{
 public:
  void rightHandSide(GridFunction& nodes, VectorGridFunction& uh, 
		     VectorGridFunction& RHS);
  // depends on position, solution to PDE, and returns RHS
  void boundaryConditions(GridFunction& nodes, VectorGridFunction& uh, 
			  VectorGridFunction& RHS)
  //depends on position and solution to PDE, 
    //returns sum of RHS and lift and external boundary conditions

};
