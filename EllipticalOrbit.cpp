#include "Orbit.h"

EllipticalOrbit::EllipticalOrbit(Coordinates coordobj){

  double rmin = params.schw.p_orb / (1.0 + params.schw.ecc);
  double rmax = params.schw.p_orb / (1.0 - params.schw.ecc);
  double xip = 0.5*(rmin+rmax);
  Sminus = params.hyperb.Sminus;
  double rstar_orb = coordobj.rstar_of_r(xip, params.schw.mass);
  
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
  p=params.schw.p_orb/params.schw.mass;
  e=params.schw.ecc;
  chi=acos(-1.0); //apastron. chi=0.0 perastron. chi=acos(0.0) start at r=params.schw.mass*p
  phi=0.0;

  //end initialize orbit
 
  cout << "p =" << p << endl;
  cout << "e =" << e << endl;
  cout << "chi =" << chi << endl;
  cout << "phi = " << phi << endl;
  cout << endl;
}

void EllipticalOrbit::dorbdchi(double &dtdchi, double &dphidchi){
  double ecoschi;
  ecoschi = e*cos(chi);
  dtdchi = pow(p, 2.0)*params.schw.mass/((p-2.0*(1.0+ecosschi))*pow((1.0+ecosschi),2.0))
    *sqrt((p-2.0*(1.0+e))*(p-2.0*(1.0-e))
	  /(p-6.0-2.0*ecoschi));
  dphidchi = sqrt(p/(p-6.0-2.0*ecoschi));
}

void EllipticalOrbit::dorbdt(){
  double dtdchi, dphidchi;
  dorbdchi(dtdchi,dphidchi);
  dchidt=1.0/dtdchi;
  dphidt=dphidchi*dchidt;
}

void EllipticalOrbit::orb_of_t(double &rp, double &drpdt, double &d2rpdt2){
  double ecosfac, mydchidt, myd2chidt2;
  ecosfac = 1.0+ e*Cos(chi);
  mydchidt = pow(ecosfac,2.0)*(-2.0+p-2.*e*cos(chi))/(params.schw.mass*(pow(p,2.0)*sqrt((-4.0*pow(e,2.0)+pow((-2.0+p),2.0)))/(-6.+p-2.*e*cos(chi))));
  myd2chidt2 = (e*pow(ecosfac,3.0)*(2.-p+2.*e*cos(chi))+(38.0+7.*pow(e,2.)+pow(p,(-19+2*p))+e*(52-11*p)*cos(chi)+7*pow(e,2.0)*cos(2*chi))*sin(chi))/(pow(params.schw.mass,2.0)*(-4.*e,2.0)+pow((-2.+p),2.)*(pow(p,4.)));
  rp=(params.schw.mass*p)/ecosfac;
  drpdt = (mydchidt*e*params.schw.mass*p*sin(chi))/pow(ecosfac,2.0);
  d2rpdt2 = (e*params.schw.mass*p*(pow(mydchidt,2.)*(2.*cos(chi)-e*(-3.+cos(2.*chi)))
		       +2.*myd2dchidt2*ecosfac*sin(chi)))
    /(2.0*pow(ecosfac,3.0));
}



