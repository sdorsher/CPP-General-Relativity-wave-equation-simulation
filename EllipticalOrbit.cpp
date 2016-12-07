#include "EllipticalOrbit.h"

EllipticalOrbit::EllipticalOrbit(){

  otype = elliptical;
 
}

void EllipticalOrbit::dorbdchi(){
  double ecoschi;
  ecoschi = e*cos(chi);
  dtdchi = pow(p, 2.0)*params.schw.mass/((p-2.0*(1.0+ecoschi))*pow((1.0+ecoschi),2.0))
    *sqrt((p-2.0*(1.0+e))*(p-2.0*(1.0-e))
	  /(p-6.0-2.0*ecoschi));
  dphidchi = sqrt(p/(p-6.0-2.0*ecoschi));
}

void EllipticalOrbit::dorbdt(){
  dorbdchi();
  dchidt=1.0/dtdchi; //in namespaces for now. should eventually become a member variable of the orbit object or this object
  dphidt=dphidchi*dchidt; // likewise
}

void EllipticalOrbit::orb_of_t(Coordinates & coords, double& rp, double& drpdt, double& d2rpdt2){
  //double rp = super.rp;
  //double drpdt = super.drpdt;
  //double d2rpdt2 = super.d2rpdt2;
  double ecosfac;
  ecosfac = 1.0+ e*cos(chi);
  dchidt = pow(ecosfac,2.0)*(-2.0+p-2.*e*cos(chi))/(params.schw.mass*(pow(p,2.0)*sqrt((-4.0*pow(e,2.0)+pow((-2.0+p),2.0)))/(-6.+p-2.*e*cos(chi))));
  d2chidt2 = (e*pow(ecosfac,3.0)*(2.-p+2.*e*cos(chi))+(38.0+7.*pow(e,2.)+pow(p,(-19+2*p))+e*(52-11*p)*cos(chi)+7*pow(e,2.0)*cos(2*chi))*sin(chi))/(pow(params.schw.mass,2.0)*(-4.*e,2.0)+pow((-2.+p),2.)*(pow(p,4.)));
  rp=(params.schw.mass*p)/ecosfac;
  drpdt = (dchidt*e*params.schw.mass*p*sin(chi))/pow(ecosfac,2.0);
  d2rpdt2 = (e*params.schw.mass*p*(pow(dchidt,2.)*(2.*cos(chi)-e*(-3.+cos(2.*chi)))
		       +2.*d2chidt2*ecosfac*sin(chi)))
    /(2.0*pow(ecosfac,3.0));
  //check here for chi value in FUTURE
}


