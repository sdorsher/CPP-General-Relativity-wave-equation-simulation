#include "Grid.h"

Grid::Grid(string fileElemBoundaries, int elemorder,int numelements):order{elemorder},NumElem{numelements},nodeLocations{0,elemorder+1,false}
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
  Array1D<double> physicalPosition(elemorder+1);
  for(int elem=0; elem<numelements; elem++)
    {
      //convert to single for loop instead of using operators on array1Ds?
      physicalPosition=((elementBoundaries[elem+1]-elementBoundaries[elem])/2.0)
        *refelem.getr()
        +((elementBoundaries[elem+1]+elementBoundaries[elem])/2.0);
      nodeLocations.append(physicalPosition);
    }

}
