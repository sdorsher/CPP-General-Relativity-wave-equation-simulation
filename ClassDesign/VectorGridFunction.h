#include "tnt_array1d.h"
#include "tnt_array2d."
#include "tnt_array1d_utils.h"
#include "tnt_array2d_utils.h"

class VectorGridFunction
{
 private:
  vector<GridFunction> data;
  int VGFvectorDim;
  int GFarrayDim;
  int GFvectorDim;

 public:
  GridFunction(int VGFvecSize, int GFvecSize, int GFarraySize, bool initZero);
  void initFromFile(string filename);
  double get(int VGFvcoord, int GFvcoord, int GFacoord);
  void set(int VGFvcoord, int GFvcoord, int GFacoord);
  int getVGFvecDim();//the dimension of the external vector
  int getGFvecDim();//the dimension of the vector within the GridFunction
  int getGFarrayDim();//the dimension of the array within the GridFunction
  VectorGridFunction(VectorGridFunction&&);
  VectorGridFunction& operator=(VectorGridFunction&&);
  //need move constructors for addition of RHS and boundary conditions
  VectorGridFunction(const VectorGridFunction&);
  VectorGridFunction& operator=(const VectorGridFunction&);
  //need copy constructors for intermediate steps of time evolution
  void operator+(VectorGridFunction, VectorGridFunction);
  //need addition operator for addition of RHS and boundary conditions
  void print(vector<string> filenames); //vector of filenames to print to
};
