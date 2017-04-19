#include "InequalityAngleConstraints.h"
#include "RigidUnitSystem.h"
#include <cmath>
#include <limits>

using namespace std;

InequalityAngleConstraints::InequalityAngleConstraints( RigidUnitSystem *rigidUnitSystem_ ) :
  rigidUnitSystem(rigidUnitSystem_)
{
}

InequalityAngleConstraints::~InequalityAngleConstraints()
{
}

double InequalityAngleConstraints::energy() {
  int nangles = angles.size();
  double Eangle = 0.0;
  double costheta;
  double diff;
  Vec3 r12;
  Vec3 r32;
  for ( int i=0; i<nangles; i++ ) {
    r12 = rigidUnitSystem->meanPositions( angles[i].p1 );
    r12 -= rigidUnitSystem->meanPositions( angles[i].p2 );
    r32 = rigidUnitSystem->meanPositions( angles[i].p3 );
    r32 -= rigidUnitSystem->meanPositions( angles[i].p2 );
    costheta = r12.dot( r32 ) / sqrt( r12.norm2()*r32.norm2() );
    
    //if angle is within limit, it contributes zero energy
    if ( costheta < angles[i].maxAllowedCosTheta ) continue;
    
    diff = acos(costheta) - angles[i].minAllowedTheta;
    Eangle += diff*diff;
  }
  return Eangle;
  
}

void InequalityAngleConstraints::addToGradient( std::vector<Vec3> &dV_dr_rigidUnitPoint,
                    std::vector<SecondDerivative> &secondDerivative_rigidUnitPoint ) {
  int nangles = angles.size();
  for ( int i=0; i<nangles; i++ ) {
    int p1 = angles[i].p1;
    int p2 = angles[i].p2;
    int p3 = angles[i].p3;
    Vec3 r12;
    Vec3 r32;
    r12 = rigidUnitSystem->meanPositions( p1 );
    r12 -= rigidUnitSystem->meanPositions( p2 );
    r32 = rigidUnitSystem->meanPositions( p3 );
    r32 -= rigidUnitSystem->meanPositions( p2 );
    double r12norm2 = r12.norm2();
    double r12norm = sqrt( r12norm2 );
    double r32norm2 = r32.norm2();
    double r32norm = sqrt( r32norm2 );
    double r12dotr32 = r12.dot( r32 );

    double costheta = r12dotr32 / ( r12norm*r32norm );
    
    //if angle is within limit, it does not contribute
    if ( costheta < angles[i].maxAllowedCosTheta ) continue;

    //compute some vector quantities we will need
    Vec3 r12perp_unitvec = r32;
    r12perp_unitvec -= r12 * (r12dotr32/r12norm2);
    r12perp_unitvec /= sqrt( r12perp_unitvec.norm2() );
    Vec3 r32perp_unitvec = r12;
    r32perp_unitvec -= r32 * (r12dotr32/r32norm2);
    r32perp_unitvec /= sqrt( r32perp_unitvec.norm2() );
    
    //compute dtheta_dr (vectors) for points r1, r2, and r3
    Vec3 dtheta_dr1 = r12perp_unitvec;
    dtheta_dr1 /= -r12norm; 
    Vec3 dtheta_dr3 = r32perp_unitvec;
    dtheta_dr3 /= -r32norm; 
    
    Vec3 dtheta_dr2 = dtheta_dr1;
    dtheta_dr2 += dtheta_dr3;
    dtheta_dr2.x = -dtheta_dr2.x;
    dtheta_dr2.y = -dtheta_dr2.y;
    dtheta_dr2.z = -dtheta_dr2.z;
    
    //now, for the three points of our angle, r1, r2, and r3,
    //compute the first and second derivatives for the corresponding
    //rigid unit points
    double prefactor = 2.0*(acos(costheta) - angles[i].minAllowedTheta);
    double prefactor_secondDerivative;
    
    const vector<int> *rup_list;
    
    rup_list = &rigidUnitSystem->getRUPlistFromP(p1);
    Vec3 dV_dr1 = dtheta_dr1;
    dV_dr1 *= prefactor/static_cast<double>(rup_list->size());
    prefactor_secondDerivative = 2.0/static_cast<double>(rup_list->size()*rup_list->size());
    for ( vector<int>::const_iterator rup_it = rup_list->begin();
          rup_it != rup_list->end();
          rup_it++ )
    {
      dV_dr_rigidUnitPoint[*rup_it] += dV_dr1;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dx2 += prefactor_secondDerivative * dtheta_dr1.x*dtheta_dr1.x;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dy2 += prefactor_secondDerivative * dtheta_dr1.y*dtheta_dr1.y;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dz2 += prefactor_secondDerivative * dtheta_dr1.z*dtheta_dr1.z;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dxdy += prefactor_secondDerivative * dtheta_dr1.x*dtheta_dr1.y;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dydz += prefactor_secondDerivative * dtheta_dr1.y*dtheta_dr1.z;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dzdx += prefactor_secondDerivative * dtheta_dr1.z*dtheta_dr1.x;
    }
    
    rup_list = &rigidUnitSystem->getRUPlistFromP(p2);
    Vec3 dV_dr2 = dtheta_dr2;
    dV_dr2 *= prefactor/static_cast<double>(rup_list->size());
    prefactor_secondDerivative = 2.0/static_cast<double>(rup_list->size()*rup_list->size());
    for ( vector<int>::const_iterator rup_it = rup_list->begin();
          rup_it != rup_list->end();
          rup_it++ )
    {
      dV_dr_rigidUnitPoint[*rup_it] += dV_dr2;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dx2 += prefactor_secondDerivative * dtheta_dr2.x*dtheta_dr2.x;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dy2 += prefactor_secondDerivative * dtheta_dr2.y*dtheta_dr2.y;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dz2 += prefactor_secondDerivative * dtheta_dr2.z*dtheta_dr2.z;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dxdy += prefactor_secondDerivative * dtheta_dr2.x*dtheta_dr2.y;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dydz += prefactor_secondDerivative * dtheta_dr2.y*dtheta_dr2.z;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dzdx += prefactor_secondDerivative * dtheta_dr2.z*dtheta_dr2.x;
    }
    
    rup_list = &rigidUnitSystem->getRUPlistFromP(p3);
    Vec3 dV_dr3 = dtheta_dr3;
    dV_dr3 *= prefactor/static_cast<double>(rup_list->size());
    prefactor_secondDerivative = 2.0/static_cast<double>(rup_list->size()*rup_list->size());
    for ( vector<int>::const_iterator rup_it = rup_list->begin();
          rup_it != rup_list->end();
          rup_it++ )
    {
      dV_dr_rigidUnitPoint[*rup_it] += dV_dr3;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dx2 += prefactor_secondDerivative * dtheta_dr3.x*dtheta_dr3.x;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dy2 += prefactor_secondDerivative * dtheta_dr3.y*dtheta_dr3.y;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dz2 += prefactor_secondDerivative * dtheta_dr3.z*dtheta_dr3.z;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dxdy += prefactor_secondDerivative * dtheta_dr3.x*dtheta_dr3.y;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dydz += prefactor_secondDerivative * dtheta_dr3.y*dtheta_dr3.z;
      secondDerivative_rigidUnitPoint[*rup_it].d2V_dzdx += prefactor_secondDerivative * dtheta_dr3.z*dtheta_dr3.x;
    }    
    
  }

}

double InequalityAngleConstraints::mismatch() {
  if ( verbose ) {
    cout << "Angles Beyond Limit:" << endl;
    cout << " point1 point2 point3 stretchAngleRadians" << endl;
  }
  int nangles = angles.size();
  int savedp1 = 0;
  int savedp2 = 0;
  int savedp3 = 0;
  double maxAngleMismatch = 0.0;
  double diff;
  double costheta;
  Vec3 r12;
  Vec3 r32;
  for ( int i=0; i<nangles; i++ ) {
    r12 = rigidUnitSystem->meanPositions( angles[i].p1 );
    r12 -= rigidUnitSystem->meanPositions( angles[i].p2 );
    r32 = rigidUnitSystem->meanPositions( angles[i].p3 );
    r32 -= rigidUnitSystem->meanPositions( angles[i].p2 );
    costheta = r12.dot( r32 ) / sqrt( r12.norm2()*r32.norm2() );
    
    //if angle is within limit, it contributes zero energy
    if ( costheta < angles[i].maxAllowedCosTheta ) continue;
    
    diff = angles[i].minAllowedTheta - acos(costheta);
    if ( verbose ) {
      cout << "  " << angles[i].p1 << " " << angles[i].p2 << " " << angles[i].p3 << " " << diff << endl;
    }
    if ( diff > maxAngleMismatch ) {
      maxAngleMismatch = diff;
      savedp1 = angles[i].p1;
      savedp2 = angles[i].p2;
      savedp3 = angles[i].p3;
    }
  }
  if ( verbose /*&& maxAngleMismatch*/ > numeric_limits<double>::epsilon() ) {
    cout << "Max Angle Beyond Limit Triple: point1 " << savedp1 << " point2 " << 
            savedp2 << " point3 " << savedp3 << " beyond-limit amount (rad) " << maxAngleMismatch << endl;
  }
  return maxAngleMismatch;
}
