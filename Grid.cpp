#include "Grid.h"

//du/dt + A du/dx + Bu = 0

Grid::Grid(int elemorder, int numelements, int nummodes):
  order{elemorder},
  NumElem{numelements},
  nodeLocs{0,elemorder + 1}, 
  refelem{elemorder},
  rschw{numelements,elemorder+1, 0.0},
  rstar{numelements,elemorder+1, 0.0},
  window{numelements,elemorder+1, 0.0},
  dwindow{numelements,elemorder+1, 0.0},
  d2window{numelements,elemorder+1, 0.0}
    
{


  //params is initialized at the beginning of this file by its 
  //inclusion in the file
  //modify such that hyperboloidal and tortoise boundaries
  // occur at element boundaries
  
  double lowerlim, upperlim; 
  
  //setup the grid and the reference element
// if (params.metric.flatspacetime) {
    // params.grid.lowerlim;
    //params.grid.upperlim;

  if(params.metric.flatspacetime){
    lowerlim=params.grid.lowerlim;
    upperlim=params.grid.upperlim;
  }
  
  if (params.metric.schwarschild) {
    lowerlim=Sminus;
    upperlim=Splus;
    cout << Sminus << " " << Splus << " set schw" << endl;
  }

  //HERE. WHY IS PARAMS>METRIC>SCHWARSZCHILD NEVER REACHED?

  //assign evenly spaced element boundaries
  for(int i = 0; i <= numelements; i++) {
    elementBoundaries.push_back(lowerlim + i * (upperlim - lowerlim) 
                                / float(numelements));
  }
  
  //Get physical positions of nodes from the reference element
  vector<double> physicalPosition(elemorder + 1);
  for(int elem = 0; elem < numelements; elem++){
    for(int node=0; node <= order; node++){
      physicalPosition[node] = ((elementBoundaries[elem + 1] 
			   - elementBoundaries[elem]) / 2.0)
	*refelem.getr()[node]
	+((elementBoundaries[elem + 1] + elementBoundaries[elem]) / 2.0);
    }
    nodeLocs.append(physicalPosition);
  }
  
  //Calculate the jacobian associated with the transformation each element
  //from the reference element to physical space
  calcjacobian();
  
}


GridFunction<double> Grid::gridNodeLocations()
{
  return nodeLocs;
}

vector<double> Grid::gridBoundaries()
{
  return elementBoundaries;
}

void Grid::calcjacobian()
{
  for(int elem = 0; elem < NumElem; elem++){
    double rx = 2.0 / (elementBoundaries[elem + 1] 
                       - elementBoundaries[elem]);
    drdx.push_back(rx);
  }
}

double Grid::jacobian(int elemnum)
{
  return drdx[elemnum];
}

int Grid::nodeOrder()
{
  return order;
}

int Grid::numberElements()
{
  return NumElem;
}

void Grid::find_extract_radii(double rfinite, double rSplus, OutputIndices& ijoutput){
  //Find the grid and element indices that correspond the the computational
  //coordinates rfinite and rSplus.
  
  bool foundfinite = false;
  bool foundSplus = false;
  for(int elem=0; elem<NumElem; elem++){
    for(int node =0; node<=order; node++){
<<<<<<< HEAD
      if((fabs(rschw.get(elem,node)-rfinite)<1.0e-5) && (!foundfinite)) { 
=======
      if((fabs(nodeLocs.get(elem,node)-rfinite)<1.0e-5) && (!foundfinite)) { 
>>>>>>> parent of 5de6444... fixed bug where location of particle was read out in different coordinate system than input. have not tested yet.
        ijoutput.ifinite = elem;
        ijoutput.jfinite = node;
        foundfinite = true;
      } else if ((fabs(nodeLocs.get(elem, node) - rSplus) < 1.0e-5) &&
                 (!foundSplus)) {
        ijoutput.iSplus = elem;
        ijoutput.jSplus = node;
        foundSplus = true;
      }


    }
   
  }
  if(!foundSplus){
    cout << "Splus output coordinate not found!"<<endl;
  }
  if(!foundfinite){
    cout << "Finite output coordinate not found!" <<endl;
  }
}
