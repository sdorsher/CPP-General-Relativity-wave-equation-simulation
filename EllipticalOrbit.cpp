#include "Orbit.h"

EllipticalOrbit::EllipticalOrbit(Coordinates coordobj){

  otype = elliptical;
  double rmin = params.schw.p_orb / (1.0 + params.schw.ecc);
  double rmax = params.schw.p_orb / (1.0 - params.schw.ecc);
  double xip = 0.5*(rmin+rmax);
  Sminus = params.hyperb.Sminus;
  double rstar_orb = coordobj.rstar_of_r(xip, params.schw.mass);
  
  double deltar = (rstar_orb - Sminus)* 2.0/params.grid.numelems; 
  coordobj.Splus = rstar_orb +round(0.5* params.grid.numelems) *deltar;
  coordobj.Rminus = rstar_orb 
    - round(0.175 * params.grid.numelems) * deltar;
  coordobj.Rplus = rstar_orb 
    + round(0.125 * params.grid.numelems) * deltar;
  coordobj.Wminus = coordobj.Rminus + params.window.noffset * deltar;
  coordobj.Wplus = coordobj.Rplus - params.window.noffset * deltar;
  
  cout << "R_star orbit" << endl;
  cout << rstar_orb << endl << endl;
  
  cout << "Sminus Rminus Rplus Splus Wminus Wplus" << endl;
  cout << coordobj.Sminus << " " << coordobj.Rminus << " " 
       << coordobj.Rplus <<" " << coordobj.Splus << " " 
       << coordobj.Wminus << " " << coordobj.Wplus << endl;
  cout << endl;

  //end initialize orbit
 
  cout << "p =" << p << endl;
  cout << "e =" << e << endl;
  cout << "chi =" << chi << endl;
  cout << "phi = " << phi << endl;
  cout << endl;
}

void EllipticalOrbit::dorbdchi(){
  double ecoschi;
  ecoschi = e*cos(chi);
  dtdchi = pow(p, 2.0)*params.schw.mass/((p-2.0*(1.0+ecosschi))*pow((1.0+ecosschi),2.0))
    *sqrt((p-2.0*(1.0+e))*(p-2.0*(1.0-e))
	  /(p-6.0-2.0*ecoschi));
  dphidchi = sqrt(p/(p-6.0-2.0*ecoschi));
}

void EllipticalOrbit::dorbdt(){
  double dtdchi, dphidchi;
  dorbdchi();
  dchidt=1.0/dtdchi; //in namespaces for now. should eventually become a member variable of the orbit object or this object
  dphidt=dphidchi*dchidt; // likewise
}

void EllipticalOrbit::orb_of_t(){
  //double rp = super.rp;
  //double drpdt = super.drpdt;
  //double d2rpdt2 = super.d2rpdt2;
  double ecosfac;
  ecosfac = 1.0+ e*Cos(chi);
  dchidt = pow(ecosfac,2.0)*(-2.0+p-2.*e*cos(chi))/(params.schw.mass*(pow(p,2.0)*sqrt((-4.0*pow(e,2.0)+pow((-2.0+p),2.0)))/(-6.+p-2.*e*cos(chi))));
  d2chidt2 = (e*pow(ecosfac,3.0)*(2.-p+2.*e*cos(chi))+(38.0+7.*pow(e,2.)+pow(p,(-19+2*p))+e*(52-11*p)*cos(chi)+7*pow(e,2.0)*cos(2*chi))*sin(chi))/(pow(params.schw.mass,2.0)*(-4.*e,2.0)+pow((-2.+p),2.)*(pow(p,4.)));
  rp=(params.schw.mass*p)/ecosfac;
  drpdt = (dchidt*e*params.schw.mass*p*sin(chi))/pow(ecosfac,2.0);
  d2rpdt2 = (e*params.schw.mass*p*(pow(dchidt,2.)*(2.*cos(chi)-e*(-3.+cos(2.*chi)))
		       +2.*d2dchidt2*ecosfac*sin(chi)))
    /(2.0*pow(ecosfac,3.0));
}



OrbitType CircularOrbit::orbType(){
  return otype;
}
