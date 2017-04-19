#include "OverlapList.h"
#include "RigidUnitSystem.h"
#include <vector>
#include "GeneralizedCoords.h"
#include <cmath>
#include <limits>
#include <iostream>

using namespace std;

OverlapList::OverlapList( 
               const RigidUnitSystem *rigidUnitSystem_,
               const VL &verletList_ ) :
    rigidUnitSystem( rigidUnitSystem_ ),
    verletList( verletList_ )
{
  update();
}

void OverlapList::update() {
  int nP = rigidUnitSystem->nPoints();
  overlapList.clear();
  #pragma omp parallel
  {
    Overlap overlap;
    int nNearbyPoints;
    const vector<VL_NeighborInfo> *nearbyPoints;
    double dist2;
    double cutoffDist;
    Vec3 delta;
    const Vec3 *p1_vec;
    const Vec3 *p2_vec;
    // variable only needed for parallel version
    #ifdef _OPENMP
    std::vector <Overlap> myOverlapList;
    #endif
    #pragma omp for
    //for ( overlap.p1 = 0; overlap.p1 < nP; overlap.p1++ ) {
    for (int p = 0; p < nP; p++) {
      overlap.p1 = p;
      nearbyPoints = &verletList.getProximityList( overlap.p1 );
      nNearbyPoints = nearbyPoints->size();
      p1_vec = &rigidUnitSystem->meanPositions(overlap.p1);
      for ( int i=0; i<nNearbyPoints; i++ ) {
        overlap.p2 = (*nearbyPoints)[i].id;
        cutoffDist = (*nearbyPoints)[i].pairCutoff;
        p2_vec = &rigidUnitSystem->meanPositions(overlap.p2);
        overlap.delta21.x = p2_vec->x - p1_vec->x;
        overlap.delta21.y = p2_vec->y - p1_vec->y;
        overlap.delta21.z = p2_vec->z - p1_vec->z;
        
        //skip this neighbor if they are too far apart.
        if ( abs(overlap.delta21.x) > cutoffDist || 
             abs(overlap.delta21.y) > cutoffDist || 
             abs(overlap.delta21.z) > cutoffDist ) continue;
        dist2 = overlap.delta21.norm2();
        if ( dist2 >= cutoffDist*cutoffDist ) continue;
        overlap.overlapDist = cutoffDist - sqrt(dist2);
    
        // if running in parallel, save a copy to private list
        #ifdef _OPENMP
        myOverlapList.push_back(overlap);
        #endif        
        //otherwise, save a copy of this overlap to shared list
        #ifndef _OPENMP 
        overlapList.push_back( overlap );
        #endif
      }
    }
    // if running in parallel, add private list to shared one
    #ifdef _OPENMP
    #pragma omp critical (OVERLAPLIST_UPDATE)
    {
      overlapList.insert(overlapList.end(),myOverlapList.begin(),myOverlapList.end());
    }
    #endif
  } // end parallel section
}

double OverlapList::energy() {  
  double Erepulsion = 0.0;  
  vector<Overlap>::const_iterator ov;
  for ( ov = overlapList.begin(); ov != overlapList.end(); ov++ ) {
    Erepulsion += ov->overlapDist * ov->overlapDist;
  }
  return Erepulsion;
}

void OverlapList::addToGradient( vector<Vec3> &dV_dr_rigidUnitPoint,
                                 vector<SecondDerivative> &secondDerivative_rigidUnitPoint ) {
  vector<Overlap>::const_iterator ov;
  for ( ov = overlapList.begin(); ov != overlapList.end(); ov++ ) {
    double deltax2 = ov->delta21.x*ov->delta21.x;
    double deltay2 = ov->delta21.y*ov->delta21.y;
    double deltaz2 = ov->delta21.z*ov->delta21.z;
    double dist2 = deltax2 + deltay2 + deltaz2;
    double dist = sqrt( dist2 );
    double cutoff_over_dist = ( ov->overlapDist + dist )/dist;
    SecondDerivative secondDerivative_unscaled;
    SecondDerivative secondDerivative;
    secondDerivative_unscaled.d2V_dx2 = 2.0 * (1.0 - cutoff_over_dist*( 1.0 - deltax2/dist2 ));
    secondDerivative_unscaled.d2V_dy2 = 2.0 * (1.0 - cutoff_over_dist*( 1.0 - deltay2/dist2 ));
    secondDerivative_unscaled.d2V_dz2 = 2.0 * (1.0 - cutoff_over_dist*( 1.0 - deltaz2/dist2 ));
    secondDerivative_unscaled.d2V_dxdy = 2.0*cutoff_over_dist/dist2*ov->delta21.x*ov->delta21.y;
    secondDerivative_unscaled.d2V_dydz = 2.0*cutoff_over_dist/dist2*ov->delta21.y*ov->delta21.z;
    secondDerivative_unscaled.d2V_dzdx = 2.0*cutoff_over_dist/dist2*ov->delta21.z*ov->delta21.x;
    Vec3 dV_dr;
    double factor1 = 2.0 * ov->overlapDist / dist;
    double factor3;

    const RigidUnitLookup::IDlist *rup_list;
    int rup;
    RigidUnitLookup::IDlist::const_iterator rup_it;
    for ( int do_twice=1; do_twice<=2; do_twice++ ) {
      if ( do_twice==1 ) {
        // First time through this loop, do the
        // gradient with respect to the degrees of freedom
        // of all rigid units corresponding to p1        
        rup_list = &rigidUnitSystem->getRUPlistFromP(ov->p1);
      }
      else {
        // For second iteration,
        // do gradient with respect to the degrees of freedom
        // of all rigid units corresponding to p2.
        factor1 = -factor1;
        rup_list = &rigidUnitSystem->getRUPlistFromP(ov->p2);
      }
      dV_dr = ov->delta21;
      dV_dr *= factor1/static_cast<double>( rup_list->size() );

      factor3 = 1.0/static_cast<double>( rup_list->size()*rup_list->size() );
      secondDerivative.d2V_dx2 = factor3*secondDerivative_unscaled.d2V_dx2;   
      secondDerivative.d2V_dy2 = factor3*secondDerivative_unscaled.d2V_dy2;   
      secondDerivative.d2V_dz2 = factor3*secondDerivative_unscaled.d2V_dz2;   
      secondDerivative.d2V_dxdy = factor3*secondDerivative_unscaled.d2V_dxdy;   
      secondDerivative.d2V_dydz = factor3*secondDerivative_unscaled.d2V_dydz;   
      secondDerivative.d2V_dzdx = factor3*secondDerivative_unscaled.d2V_dzdx;
      for ( rup_it = rup_list->begin();
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


double OverlapList::mismatch() {
  double maxOverlap = 0.0;
  if ( verbose ) {
    cout << "Overlaps:" << endl;
    cout << " point1 point2 overlapDist" << endl;
  }
  int p1 = 0;
  int p2 = 0;
  vector<Overlap>::const_iterator ov;
  for ( ov = overlapList.begin(); ov != overlapList.end(); ov++ ) {
    if ( verbose ) {
      cout << "  " << ov->p1 << " " << ov->p2 << " " << ov->overlapDist << endl;
    }
    if ( ov->overlapDist > maxOverlap ) {
      maxOverlap = ov->overlapDist;
      p1 = ov->p1;
      p2 = ov->p2;
    }
  }
  if ( verbose && maxOverlap > numeric_limits<double>::epsilon() ) {
    double sep = sqrt((rigidUnitSystem->meanPositions(p1)-rigidUnitSystem->meanPositions(p2)).norm2());
    cout << "Max Overlap: point1 " << p1 << " point2 " << p2 << 
    " separation " << sep <<
    " overlapDist " << maxOverlap << 
    " constraint " << sep + maxOverlap << endl;
  }
  return maxOverlap;
}

