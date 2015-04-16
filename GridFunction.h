#ifndef _GRIDFUNCTION_H
#define _GRIDFUNCTION_H
#include "tnt/tnt.h"
#include "TNT2.h"
#include <fstream>

using namespace std;

class GridFunction
{
 public:
  GridFunction(int vecSize, int arraySize,bool initZeros);
  void set(int vcoord, TNT::Array1D<double>);
  void set(int vcoord, int acoord,double value);
  void append(TNT::Array1D<double> array);
  double get(int vcoord, int acoord);
  TNT::Array1D<double> get(int vcoord);
  int gridDim();
  int pointsDim();
  void save(string filename);
  
  // private:
 private:
  vector<TNT::Array1D<double>> data;
  int GFvectorDim;
  int GFarrayDim;

  

//  void initFromFile(string filename);






//  GridFunction(GridFunction&&);
//  GridFunction& operator=(GridFunction&&);
  //need move constructors for addition of RHS and boundary conditions
//  GridFunction(const GridFunction&);
//  GridFunction& operator=(const GridFunction&);
  //need copy constructors for intermediate steps of time evolution
 
};

GridFunction operator+(GridFunction gf1,GridFunction gf2);
GridFunction operator*(double A,GridFunction gf2);


#endif
