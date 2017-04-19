#include "Ropes.h"
#include "RigidUnitSystem.h"
#include <vector>
#include "GeneralizedCoords.h"
#include <cmath>
#include <limits>
#include <iostream>


using namespace std;

Ropes::Ropes( const RopesAssignment& ropes_,
              const RigidUnitSystem *rigidUnitSystem_ ) :
  ropes( ropes_ ),
  rigidUnitSystem( rigidUnitSystem_ )
{
}


double Ropes::energy() {
  int nropes = ropes.size();
  double Erope = 0.0;
  #pragma omp parallel reduction(+ : Erope)
  {   
    int p1;
    int p2;
    double ropelength;
    double dist2;
    double diff;
    #pragma omp for
    for ( int i=0; i<nropes; i++ ) {
      p1 = ropes[i].p1;
      p2 = ropes[i].p2;
      ropelength = ropes[i].length;
      dist2 = rigidUnitSystem->meanPositions(p1).dist2( rigidUnitSystem->meanPositions(p2) );
      if ( dist2 < ropelength*ropelength ) continue;
      diff = sqrt(dist2) - ropelength;
      //cout << p1 << " " << p2 << " " << sqrt(dist2) << endl;
      Erope += diff*diff;
    }
  } // end parallel section
  return Erope;
}

void Ropes::addToGradient( vector<Vec3> &dV_dr_rigidUnitPoint,
                           vector<SecondDerivative> &secondDerivative_rigidUnitPoint ) {
  int nropes = ropes.size();
  for ( int i=0; i<nropes; i++ ) {
    int p1 = ropes[i].p1;
    int p2 = ropes[i].p2;
    double ropelength = ropes[i].length;

    Vec3 delta = rigidUnitSystem->meanPositions(p2);
    delta -= rigidUnitSystem->meanPositions(p1);

    double deltax2 = delta.x*delta.x;
    double deltay2 = delta.y*delta.y;
    double deltaz2 = delta.z*delta.z;
    double dist2 = deltax2 + deltay2 + deltaz2;
    if ( dist2 < ropelength*ropelength ) continue;
    // otherwise, dist2 is greater than the ropelength2
    // (the rope is stretched and must pull the points towards each other)

    double dist = sqrt( dist2 );
    double cutoff_over_dist = ropelength/dist;
    SecondDerivative secondDerivative_unscaled;
    SecondDerivative secondDerivative;
    secondDerivative_unscaled.d2V_dx2 = 2.0 * (1.0 - cutoff_over_dist*( 1.0 - deltax2/dist2 ));
    secondDerivative_unscaled.d2V_dy2 = 2.0 * (1.0 - cutoff_over_dist*( 1.0 - deltay2/dist2 ));
    secondDerivative_unscaled.d2V_dz2 = 2.0 * (1.0 - cutoff_over_dist*( 1.0 - deltaz2/dist2 ));
    secondDerivative_unscaled.d2V_dxdy = 2.0*cutoff_over_dist/dist2*delta.x*delta.y;
    secondDerivative_unscaled.d2V_dydz = 2.0*cutoff_over_dist/dist2*delta.y*delta.z;
    secondDerivative_unscaled.d2V_dzdx = 2.0*cutoff_over_dist/dist2*delta.z*delta.x;

    Vec3 dV_dr;
    double factor1 = 2.0 * (ropelength - dist) / dist;
    double factor3;
    
    for ( int do_twice=1; do_twice<=2; do_twice++ ) {
      const vector<int> *rup_list;
      if ( do_twice==1 ) {
        // First time through this loop, do the
        // gradient with respect to the degrees of freedom
        // of all rigid units corresponding to p1
        rup_list = &rigidUnitSystem->getRUPlistFromP(p1);
      }
      else {
        // For second iteration,
        // do gradient with respect to the degrees of freedom
        // of all rigid units corresponding to p2.
        factor1 = -factor1;
        rup_list = &rigidUnitSystem->getRUPlistFromP(p2);
      }
      dV_dr = delta;
      dV_dr *= factor1/static_cast<double>( rup_list->size() );

      factor3 = 1.0/static_cast<double>( rup_list->size()*rup_list->size() );
      secondDerivative.d2V_dx2 = factor3*secondDerivative_unscaled.d2V_dx2;   
      secondDerivative.d2V_dy2 = factor3*secondDerivative_unscaled.d2V_dy2;   
      secondDerivative.d2V_dz2 = factor3*secondDerivative_unscaled.d2V_dz2;   
      secondDerivative.d2V_dxdy = factor3*secondDerivative_unscaled.d2V_dxdy;   
      secondDerivative.d2V_dydz = factor3*secondDerivative_unscaled.d2V_dydz;   
      secondDerivative.d2V_dzdx = factor3*secondDerivative_unscaled.d2V_dzdx;
      int rup;
      for ( vector<int>::const_iterator rup_it = rup_list->begin();
            rup_it != rup_list->end();
            rup_it++ )
      {
        rup = *rup_it;
        dV_dr_rigidUnitPoint[rup] += dV_dr;
  
        secondDerivative_rigidUnitPoint[rup].d2V_dx2 += secondDerivative.d2V_dx2;
        secondDerivative_rigidUnitPoint[rup].d2V_dy2 += secondDerivative.d2V_dy2;
        secondDerivative_rigidUnitPoint[rup].d2V_dz2 += secondDerivative.d2V_dz2;
        secondDerivative_rigidUnitPoint[rup].d2V_dxdy += secondDerivative.d2V_dxdy;
        secondDerivative_rigidUnitPoint[rup].d2V_dydz += secondDerivative.d2V_dydz;
        secondDerivative_rigidUnitPoint[rup].d2V_dxdy += secondDerivative.d2V_dxdy;
      }
    }
  }
}

double Ropes::mismatch() {
  if ( verbose ) {
    cout << "Ropes Stretched:" << endl;
    cout << " point1 point2 stretchDist" << endl;
  }
  int nropes = ropes.size();
  int savedp1 = 0;
  int savedp2 = 0;
  double maxRopeMismatch = 0.0;
  int p1;
  int p2;
  double ropelength;
  double dist2;
  double diff;
  for ( int i=0; i<nropes; i++ ) {
    p1 = ropes[i].p1;
    p2 = ropes[i].p2;
    ropelength = ropes[i].length;
    dist2 = rigidUnitSystem->meanPositions(p1).dist2( rigidUnitSystem->meanPositions(p2) );
    if ( dist2 < ropelength*ropelength ) continue;
    diff = sqrt(dist2) - ropelength;
    if ( verbose ) {
      cout << "  " << p1 << " " << p2 << " " << diff << endl;
    }
    if ( diff > maxRopeMismatch ) {
      maxRopeMismatch = diff;
      savedp1 = p1;
      savedp2 = p2;
    }
  }
  if ( verbose && maxRopeMismatch > numeric_limits<double>::epsilon() ) {
    cout << "Max Rope Over-stretched Pair: point1 " << savedp1 << " point2 " << 
            savedp2 << " over-stretch amount " << maxRopeMismatch << endl;
  }
  return maxRopeMismatch;
}

