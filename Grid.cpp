#include "Grid.h"

//du/dt + A du/dx + Bu = 0

Grid::Grid(int elemorder, int numelements, int nummodes, double lowerlim, 
           double upperlim):
  order{elemorder},
  NumElem{numelements},
  nodeLocs{0,elemorder + 1}, 
  refelem{elemorder},
  rschw{numelements,elemorder+1},
  rstar{numelements,elemorder+1}
{


  //params is initialized at the beginning of this file by its 
  //inclusion in the file
  //modify such that hyperboloidal and tortoise boundaries
  // occur at element boundaries

  //setup the modes
  Modes lmmodes(params.modes.lmax);

  if(params.metric.schwarschild){
    double rmin = params.schw.p_orb / (1.0 + params.schw.ecc);
    double rmax = params.schw.p_orb / (1.0 - params.schw.ecc);
    double xip = 0.5*(rmin+rmax);
    Sminus = params.hyperb.Sminus;
    double rstar_orb = rstar_of_r(xip, params.schw.mass);
    
    double deltar = (rstar_orb - Sminus)* 2.0/params.grid.numelems; 
    Splus = rstar_orb +round(0.5* params.grid.numelems) *deltar;
    Rminus = rstar_orb 
      - round(0.175 * params.grid.numelems) * deltar;
    Rplus = rstar_orb 
      + round(0.125 * params.grid.numelems) * deltar;
    Wminus = Rminus + params.window.noffset * deltar;
    Wplus = Rplus - params.window.noffset * deltar;

    cout << "R_star orbit" << endl;
    cout << rstar_orb << endl << endl;
    
    cout << "Sminus Rminus Rplus Splus Wminus Wplus" << endl;
    cout << Sminus << " " << Rminus << " " 
	 << Rplus <<" " << Splus << " " 
	 << Wminus << " " << Wplus << endl;
    cout << endl;


    //initialize orbit
    initialize_orbit();
    cout << "p =" << p << endl;
    cout << "e =" << e << endl;
    cout << "chi =" << chi << endl;
    cout << "phi = " << phi << endl;
    cout << endl;

 
    if(params.opts.useSource) {
      init_source( lmmodes, params.schw.mass);
    }
    R1= invert_tortoise(Rminus, params.schw.mass) + 2.0*params.schw.mass;
    R2 = invert_tortoise(Rplus, params.schw.mass) + 2.0* params.schw.mass;
    w1 = params.schw.p_orb-(invert_tortoise(2.0*deltar, params.schw.mass)
                          +2.0*params.schw.mass)-R1;
    w2 = R2 - (params.schw.p_orb + invert_tortoise(2.0*deltar, params.schw.mass)+2.0*params.schw.mass);
    nmodes = lmmodes.ntotal;


    cout << "R1 R2 w1 w2" << endl;
    cout << R1 << " " << R2 << " " << w1 <<  " " << w2 << endl << endl;

    
    if(params.opts.useSource) {
      set_window(R1, w1, 1.0, 1.5, R2, w2, 1.0, 1.5, lmmodes.ntotal);
    }
  }
  
  double lowlim, uplim; 
  
  //setup the grid and the reference element
  if (params.metric.flatspacetime) {
    lowlim = params.grid.lowerlim;
    uplim = params.grid.upperlim;
    
   } else if (params.metric.schwarschild) {
    lowlim = Sminus;
    uplim = Splus;
  }



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
      if((fabs(nodeLocs.get(elem,node)-rfinite)<1.0e-5) && (!foundfinite)) { 
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
}
