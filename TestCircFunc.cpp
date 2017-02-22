#include "Orbit.h"
#include "CircularOrbit.h"
#include "globals.h"
#include "ConfigParams.h"

void orbitfn(Orbit* orb);

int main(){
  
 CircularOrbit* corb = new CircularOrbit();
 cout << corb->p << " " << corb->e << " " << corb->chi << " " << corb->phi << " " << corb->E << " " << corb->L << endl;

 corb->phi=corb->phi_of_t(0.1);
 cout <<  corb->circ_E() << endl;
 cout << corb->circ_L() << endl;
 cout << corb->p << " " << corb->e << " " << corb->chi << " " << corb->phi << " " << corb->E << " " << corb->L << endl;

 orbitfn(corb);
 

 delete corb;
}

void orbitfn(Orbit* orb){
  if((orb->orbType())==circular){
    CircularOrbit* corb2=dynamic_cast<CircularOrbit*>(orb);
    corb2->phi_of_t(200);
    cout << corb2->phi << endl;
    cout << corb2->circ_E() << endl;
    cout << corb2->circ_L() << endl;
    cout << orb->p << " " << orb->e << " " << orb->chi << " " << orb->phi << " " << orb->E << " " << orb->L << endl;
  }
}
