#include "tnt/tnt.h"
#include "TNT2.h"

using namespace std;

class GridFunction
{
 public:
  GridFunction(int vecSize, int arraySize,bool initZeros);
  void set(int vcoord, TNT::Array1D<double>);
  void set(int vcoord, int acoord,double value);
  void append(TNT::Array1D<double> array);
  double get(int vcoord, int acoord);
  int gridDim();
  int functionDim();

  //  void operator+(GridFunction, GridFunction);
  //need addition operator for addition of RHS and boundary conditions


  // private:
 public:
  vector<TNT::Array1D<double>> data;
  int GFvectorDim;
  int GFarrayDim;
  /* void initFromFile(string filename);





  void save(string filename)


//  GridFunction(GridFunction&&);
//  GridFunction& operator=(GridFunction&&);
  //need move constructors for addition of RHS and boundary conditions
//  GridFunction(const GridFunction&);
//  GridFunction& operator=(const GridFunction&);
  //need copy constructors for intermediate steps of time evolution
  */
};

