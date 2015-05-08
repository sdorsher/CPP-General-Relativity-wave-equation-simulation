#include "Grid.h"

Grid::Grid(string fileElemBoundaries, int elemord,int numelements):order{elemord},NumElem{numelements},nodeLocs{0,elemord+1,false}
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

  ReferenceElement refelem(elemord);
  //when we generalize this, use map to store order, element pairs
  //so they do not need to be recalculated with each element of the same
  //order
  
  Array1D<double> physicalPosition(elemord+1);
  for(int elem=0; elem<numelements; elem++)
    {
      physicalPosition=((elementBoundaries[elem+1]-elementBoundaries[elem])/2.0)
        *refelem.getr()
        +((elementBoundaries[elem+1]+elementBoundaries[elem])/2.0);
      nodeLocs.append(physicalPosition);
    }

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
      double rx= (elementBoundaries[elem+1]-elementBoundaries[elem])/2.0;
      //cout << elem << " " <<rx << endl;
      drdx.push_back(rx);
    }
}

vector<double> Grid::jacobian()
{
  return drdx;
}
