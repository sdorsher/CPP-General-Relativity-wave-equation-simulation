#include "Modes.h"

Modes::Modes(int lmax) 
{
  set_lm_mode_info(lmax);
}

int Modes::n_of_l(int l)
{
  if (l % 2 == 0) {
    return (l + 2) / 2;
  } else {
    return (l + 1) / 2;
  }
}

int Modes::nmodes_of_l(int lmax)
{
  int nmodes = 0;
  
  for(int l=0; l<=lmax; l++) {
    nmodes += n_of_l(l);
  }
  return nmodes;
}

void Modes::set_lm_mode_info(int lmax) {
  ntotal = nmodes_of_l(lmax);
  ll.resize(ntotal);
  mm.resize(ntotal);
  psil.resize(lmax+1);
  psitl.resize(lmax+1);
  psiphil.resize(lmax+1);
  psirl.resize(lmax+1);
  int n = 0;
  for(int l=0; l<=lmax; l++){
    for(int m = l%2; m<=l; m+=2){
      ll[n]=l;
      mm[n]=m;
      n++;
    }
  }
}

void Modes::sum_m_modes(TwoDVectorGridFunction<complex<double>> uh,double time,int index1,int index2, Orbit * orb)
{
  double m_fold_factor, y_lm, omega, ecosfac, radius;
  complex<double> phase, phifactor;
  complex<double> zi{0.,1.0};


  for(int i=0; i<psil.size(); i++){
    psil.at(i)= 0.0;
    psitl.at(i)=0.0;
    psiphil.at(i)=0.0;
    psirl.at(i)=0.0;
  }
  ecosfac = 1.0 + orb->e*cos(orb->chi);
  radius = (params.schw.mass*orb->p)/ecosfac;
  //if(orb->orbType()==circular){
  // CircularOrbit * corb = dynamic_cast<CircularOrbit *>(orb);
  // orb->phi=corb->phi_of_t(time);
  //}
  for (int i=0; i<mm.size(); i++){
    if(mm.at(i)==0) {
      m_fold_factor = 1.0;
    } else {
      m_fold_factor = 2.0;
    }
    phifactor = zi*mm.at(i);
    complex<double> tempphase(cos(mm.at(i)*orb->phi), sin(mm.at(i)*orb->phi));
    phase = tempphase;
    y_lm = gsl_sf_legendre_sphPlm(ll.at(i), mm.at(i), 0.0);
    //y_lm = legendre_sphPlm(ll.at(i), mm.at(i), 0.0);
    psil.at(ll.at(i))= psil.at(ll.at(i))
      +m_fold_factor * y_lm * (phase.real() *uh.get(i,0,index1,index2).real());
    psitl.at(ll.at(i))=psitl.at(ll.at(i))
      +m_fold_factor*y_lm* (phase.real()* uh.get(i,2,index1,index2).real());
    psiphil.at(ll.at(i)) = psiphil.at(ll.at(i))
      +m_fold_factor * y_lm* (-phifactor.imag()*phase.real()*uh.get(i,0,index1,index2).imag());
    psirl.at(ll.at(i)) = psirl.at(ll.at(i))
      +m_fold_factor * y_lm *(phase.real()*uh.get(i, 1, index1, index2).real());
  }

  //convert to physical modes. the time and radial derivatives have not been checked in Fortran
  for(int i = 0; i < psil.size(); i++){
    psirl.at(i) = psirl.at(i)/(radius - 2.0* params.schw.mass)
      - psil.at(i)/pow(radius, 2.0);
    psil.at(i) = psil.at(i)/radius;
    psitl.at(i) = psitl.at(i) /radius;
    psiphil.at(i)= psiphil.at(i)/radius;
  }

}
