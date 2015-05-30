#include "DiffEq.h"

GridFunction<Array2D<double>> setupAmatrix(GridFunction<double>& nodes)
{
  GridFunction<Array2D<double>> gf(nodes.gridDim(),nodes.pointsDim());
  for(int i=0; i<nodes.gridDim(); i++){
    for(int j=0; j<nodes.pointsDim(); j++){
      Array2D<double> A(3,3,0.0);
      A[1][2]=-pow(params.waveeq.speed,2.0);
      A[2][1]=-1.0;
      gf.set(i,j,A);
    }
  }
  return gf;
}

GridFunction<Array2D<double>> setupBmatrix(GridFunction<double>& nodes)
{
  GridFunction<Array2D<double>> gf(nodes.gridDim(),nodes.pointsDim());
  for(int i=0; i<nodes.gridDim(); i++){
    for(int j=0; j<nodes.pointsDim(); j++){
      Array2D<double> B(3,3,0.0);
      B[0][1]=1.0;
      gf.set(i,j,B);
    }
  }
  return gf;
}
