#include "SharedPoints.h"
#include "RigidUnitSystem.h"
#include <vector>
#include "GeneralizedCoords.h"
#include <cmath>
#include <limits>
#include <iostream>

using namespace std;

SharedPoints::SharedPoints( const RigidUnitSystem *rigidUnitSystem_ ) :
  rigidUnitSystem( rigidUnitSystem_ )
{
  int nP = rigidUnitSystem->nPoints();
  for ( int p = 0; p < nP; p++ ) {
    if ( rigidUnitSystem->getRUPlistFromP(p).size() > 1 ) {
      sharedPoints.push_back( p );
    }
  }
  std::vector<int>(sharedPoints).swap(sharedPoints);
}

SharedPoints::~SharedPoints()
{
}

double SharedPoints::energy() {
  double EsharedPoint = 0.0;
  #pragma omp parallel for reduction(+ : EsharedPoint)
  for ( size_t i = 0; i < sharedPoints.size(); i++ ) {
    size_t p = sharedPoints[i];
    const RigidUnitLookup::IDlist *rup_list = 
      &rigidUnitSystem->getRUPlistFromP(p);
    double nrup_double = static_cast<double>( rup_list->size() );
       const Vec3 *meanvec = &rigidUnitSystem->meanPositions(p);
    for ( vector<int>::const_iterator rup_it = rup_list->begin();
          rup_it != rup_list->end();
          rup_it++ )
    {
      EsharedPoint += nrup_double * meanvec->dist2( rigidUnitSystem->absolutePositions(*rup_it) );
    }

  }
  return EsharedPoint;
}

void SharedPoints::addToGradient( vector<Vec3> &dV_dr_rigidUnitPoint,
                                  vector<SecondDerivative> &secondDerivative_rigidUnitPoint ) {
  #pragma omp parallel for
  for ( size_t i = 0; i < sharedPoints.size(); i++ ) {
    size_t p = sharedPoints[i];
    const Vec3 *meanvec = &rigidUnitSystem->meanPositions(p);
    const vector<int> *rup_list = 
      &rigidUnitSystem->getRUPlistFromP(p);
    double nrup_double = static_cast<double>( rup_list->size() );
    double factor1 = 2.0*nrup_double;
    double secondDerivative = factor1 - 2.0; // same as 2*( nrup - 1 )
    int rup;
    // No need to use a critical directive here when parallelizing 
    // (the way I did with Ropes and OverlapList)
    // because each RUP belongs to only one point, so those variables will only ever
    // be written by a single iteration of the loop.
    for ( vector<int>::const_iterator rup_it = rup_list->begin();
          rup_it != rup_list->end();
          rup_it++ )
    {
      rup = *rup_it;
      Vec3 delta = rigidUnitSystem->absolutePositions(rup);
      delta -= *meanvec;
      delta *= factor1;
      
      dV_dr_rigidUnitPoint[rup] += delta;
      secondDerivative_rigidUnitPoint[rup].d2V_dx2 += secondDerivative;
      secondDerivative_rigidUnitPoint[rup].d2V_dy2 += secondDerivative;
      secondDerivative_rigidUnitPoint[rup].d2V_dz2 += secondDerivative;
    }
  }
}

double SharedPoints::mismatch() {
  double maxMismatch2 = 0.0;
  size_t savedp = 0;
  double dist2;
  for ( size_t i = 0; i < sharedPoints.size(); i++ ) {
    size_t p = sharedPoints[i];
    const vector<int> *rup_list = 
      &rigidUnitSystem->getRUPlistFromP(p);
    const Vec3 *meanvec = &rigidUnitSystem->meanPositions(p); 
    for ( vector<int>::const_iterator rup_it = rup_list->begin();
          rup_it != rup_list->end();
          rup_it++ )
    {
      dist2 = meanvec->dist2( rigidUnitSystem->absolutePositions(*rup_it) );
      if ( dist2 > maxMismatch2 ) {
        maxMismatch2 = dist2;
        savedp = p;
	}
    }
  } // end loop
  if ( verbose && maxMismatch2 > numeric_limits<double>::epsilon() ) {
    cout << "Max Shared Point to Mean Point: point " << 
            savedp << " dist " << sqrt( maxMismatch2 ) << endl;
  }

  return sqrt( maxMismatch2 );
}

