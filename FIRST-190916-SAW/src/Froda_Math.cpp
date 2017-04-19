#include "global_defs.h"
#include "generalUtils.h"
#include "interact.h"
#include "Froda.h"
#include "mt19937ar.h"
#include "flexweb.h"

#include <signal.h>

extern Parameters parameters;


////////////////////////////////////////////////////////////////////////////////
// Description: generates a unit vector in a random direction.
////////////////////////////////////////////////////////////////////////////////
Vector Froda::randomUnitVector( ){
  double myphi, mycostheta, mysintheta, range;
  Vector res(0,0,0);
      range = TWOPI;
      myphi = range*genrand_real1(); // 0 to 2pi
      range = 2;
      mycostheta = (range*genrand_real1()) - (double) 1.0; //-1 to 1
      mysintheta = sqrt(1 - pow(mycostheta,2));

      res.x = mysintheta * cos(myphi);
      res.y = mysintheta * sin(myphi); 
      res.z = mycostheta;
 
  return res;
} 



////////////////////////////////////////////////////////////////////////////////
// Description: guess a rotor for a ghost
// Not currently in use, but might become useful; keep the algebra.
////////////////////////////////////////////////////////////////////////////////
Vector Froda::estimateRotor( int ghostid){
    int ghostsize = ghost.at(ghostid).bondVector.size();
    double A0, A1, A2, A3, A4, A5, A6, A7, A8;
    A0 = A1 = A2 = A3 = A4 = A5 = A6 = A7 = A8 = 0;
    for ( int p = 0; p < ghostsize; p++){
      int theatom = ghost[ghostid].atoms.at(p);
      double aweight = frodaAtom.at(theatom).weight;
      A0 += ( aweight * ghost[ghostid].posAtom.at(p).x * ghost[ghostid].bondVector.at(p).x );
      A1 += ( aweight * ghost[ghostid].posAtom.at(p).x * ghost[ghostid].bondVector.at(p).y );
      A2 += ( aweight * ghost[ghostid].posAtom.at(p).x * ghost[ghostid].bondVector.at(p).z );
      A3 += ( aweight * ghost[ghostid].posAtom.at(p).y * ghost[ghostid].bondVector.at(p).x );
      A4 += ( aweight * ghost[ghostid].posAtom.at(p).y * ghost[ghostid].bondVector.at(p).y );
      A5 += ( aweight * ghost[ghostid].posAtom.at(p).y * ghost[ghostid].bondVector.at(p).z );
      A6 += ( aweight * ghost[ghostid].posAtom.at(p).z * ghost[ghostid].bondVector.at(p).x );
      A7 += ( aweight * ghost[ghostid].posAtom.at(p).z * ghost[ghostid].bondVector.at(p).y );
      A8 += ( aweight * ghost[ghostid].posAtom.at(p).z * ghost[ghostid].bondVector.at(p).z );
    }
    double C1, C2, C3, C4, C5, C6, C7, C8, C9;
    C1 = A0 + A4;
    C2 = A4 + A8;
    C3 = A8 + A0;

    C4 = A7 - A5;
    C5 = A2 - A6;
    C6 = A3 - A1;

    C7 = 0.5 * ( A1 + A3 ); 
    C8 = 0.5 * ( A5 + A7 ); 
    C9 = 0.5 * ( A6 + A2 ); 

    double D = (C1 * C7 * C7)  + (C2 * C8 * C8) + (C3 * C9 * C9) +( 2*C7 * C8 * C9) - (C1 * C2 * C3); 
 
    // <- get the estimated rotor here
    Vector rotor_estimate;

    rotor_estimate.x = ( 1.0/D) * ( ( - C1*C3 + C8*C8)*C4   - (C1*C7 + C8*C9)*C5  - (C3*C9 + C7*C8)*C6 );
    rotor_estimate.y = ( 1.0/D) * ( ( - C2*C1 + C9*C9)*C5   - (C2*C8 + C9*C7)*C6  - ( C1*C7 +C8*C9)*C4 );
    rotor_estimate.z = ( 1.0/D) * ( ( - C3*C2 + C7*C7)*C6   - ( C3*C9 + C7*C8)*C4  - (C2*C8 + C9*C7)*C5 );
  return rotor_estimate;
}



