#include "tnt_array1d.h"
#include "tnt_array2d."
#include "tnt_array1d_utils.h"
#include "tnt_array2d_utils.h"

class GridFunction
{
 private:
  vector<Array1D<double>> data;
  int GFvectorDim;
  int GFarrayDim;
 public:
  GridFunction(int vecSize, int arraySize,bool initZeros);
  void initFromFile(string filename);
  double get(int vcoord, int acoord);
  double set(int vcoord, int acoord);
  int getGFvecDim();
  int getGFarrDim();  
  GridFunction(GridFunction&&);
  GridFunction& operator=(GridFunction&&);
  //need move constructors for addition of RHS and boundary conditions
  GridFunction(const GridFunction&);
  GridFunction& operator=(const GridFunction&);
  //need copy constructors for intermediate steps of time evolution
  void operator+(GridFunction, GridFunction);
  //need addition operator for addition of RHS and boundary conditions
  void print(string filename)
};

