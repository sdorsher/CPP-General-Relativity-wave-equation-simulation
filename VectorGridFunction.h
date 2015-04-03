#include "tnt.h"
#include "TNT2.h"
#include "GridFunction.h"
#include <fstream>


class VectorGridFunction
{
 public:
  vector<GridFunction> data;
  int VGFvectorDim;
  int GFarrayDim;
  int GFvectorDim;

 public:
  VectorGridFunction(int VGFvecSize, int GFvecSize, int GFarraySize, bool initZero);
  double get(int VGFvcoord, int GFvcoord, int GFacoord);
  Array1D<double> get(int VGFvcoord,int GFvcoord);
  GridFunction get(int VGFvcoord);
  void set(int VGFvcoord, int GFvcoord, int GFacoord,double value);
  void set(int VGFcoord,GridFunction value);
  void set(int VGFcoord,int GFcoord,TNT::Array1D<double> arr);
  int vectorDim();//the dimension of the external vector
  int gridDim();//the dimension of the vector within the GridFunction
  int pointsDim();//the dimension of the array within the GridFunction
  void append(GridFunction gf);
  void save(string filenames); //vector of filenames to print to
  //  void initFromFile(string filename);
  vector<double> getVector(int GFvcoord, int GFacoord);
  
};

VectorGridFunction operator+(VectorGridFunction, VectorGridFunction);
  //need addition operator for addition of RHS and boundary conditions

using FUNCTYPE = void(GridFunction&, VectorGridFunction&, VectorGridFunction&, int k, int i);


void loop(GridFunction& grid, VectorGridFunction& uh, VectorGridFunction& RHS, FUNCTYPE func);

void testfunc(GridFunction& grid, VectorGridFunction& uh, VectorGridFunction& RHS, int k, int i);
