#include "Grid.h"

Grid::Grid(string fileElemBoundaries, int elemorder,int numelements):order{elemorder},NumElem{numelements},nodeLocs{0,elemorder+1,false},drdx{NULL}
{
  ifstream fs;
  fs.open(fileElemBoundaries);
  double data;
  fs >> data;
  while(!fs.eof())
    {
      elementBoundaries.push_back(data);
      fs >> data;
    }
  fs.close();
  if(elementBoundaries.size()-1 != numelements)
    {
      throw invalid_argument("Element boundaries too long or too short for number of elements given.");
    }

  ReferenceElement refelem(elemorder);
  //when we generalize this, use map to store order, element pairs
  //so they do not need to be recalculated with each element of the same
  //order
  

  calcjacobian();
  

}

int Grid::numberElements()
{
  return NumElem;
}

GridFunction Grid::gridNodeLocations()
{
  return nodeLocs;
}

vector<double> Grid::gridBoundaries()
{
  return elementBoundaries;
}

void Grid::calcjacobian()
{
  for(int elem=0; elem<NumElem; elem++)
    {

      drdx.push_back((elementBoundaries[elem+1]-elementBoundaries[elem])/2.0);
    }
}

vector<double> Grid::jacobian()
{
  return drdx;
}
