#include "tnt_array1d.h"
#include "tnt_array2d."
#include "tnt_array1d_utils.h"
#include "tnt_array2d_utils.h"
#include "ReferenceElement.h"
#include "GridFunction.h"

class Grid
{

 private:
  string fileElemBoundaries;
  int NumElem;
  Array1D elementBoundaries; //2N
  GridFunction nodesLocations; //NxNp

 public:
  Grid(string fileElemBoundaries, ReferenceElement refelem);
  //initializes elementBoundaries from file,
  //obtains node locations from reference element and puts them in array

  int getN();//return number of elements, calculated from input file
  GridFunction getNodeLocations();  
};


