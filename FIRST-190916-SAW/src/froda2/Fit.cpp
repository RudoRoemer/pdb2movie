#include "Fit.h"
#include "Rotator.h"
#include "MultiVarFunction_Adapter_Fit.h"
#include <cmath>
#include <limits>

using namespace std;

Fit::Fit() :
  func( NULL ),
  rotormagcutoff2( 1.0*1.0 )
{
  func = new MultiVarFunction_Adapter_Fit( this );
  targetAbsoluteCenter = Vec3(0,0,0);

  minimizer.setLimitingParabolicMagnification( 10.0 );
  minimizer.set_requireRigorousLineMinimization( false );
  minimizer.setTrialStepDQmax( 0.1 );
  minimizer.setTrialStepLambda( 0.01 );
  rotorIsZero = true;
  hasZeroRadius = false;
  numberOfSavedRotations = 0;
}

Fit::~Fit()
{
}

void Fit::setSourceBasePoints( const std::vector<Vec3> &base_, double radiusGyration_ ) {
  base = base_;
  rotor = Vec3( 0,0,0 );
  rotated = base;
  radiusGyration = radiusGyration_;

  rotorIsZero = true;  
  hasZeroRadius = radiusGyration <= numeric_limits<double>::epsilon();
  isUpdatedE = false;
  isUpdatedGrad = false;
  isUpdated2ndDeriv = false;
  numberOfSavedRotations = 0;
}

void Fit::setSourceAbsolutePoints( const std::vector<Vec3> &absolute_ ) {
  Vec3 center(0,0,0); 
  double rg = 0;
  size_t N = absolute_.size();
  for ( size_t i=0; i < N; i++ ) { 
    center += absolute_[i];
  }
  center /= N;

  base.resize( N );
  for ( size_t i=0; i < N; i++ ) { 
    base[i] = absolute_[i] - center;
    rg += base[i].norm2();
  }  
  rg = sqrt( rg/N );
  
  rotor = Vec3( 0,0,0 );
  rotated = base;
  radiusGyration = rg;

  rotorIsZero = true;
  hasZeroRadius = radiusGyration <= numeric_limits<double>::epsilon();
  isUpdatedE = false;
  isUpdatedGrad = false;
  isUpdated2ndDeriv = false;
  numberOfSavedRotations = 0;

}

void Fit::setTargetAbsolutePoints( const std::vector<Vec3> &absolute_ ) {
  targetAbsoluteCenter = Vec3(0,0,0); 
  size_t N = absolute_.size();
  for ( size_t i=0; i < N; i++ ) { 
    targetAbsoluteCenter += absolute_[i];
  }
  targetAbsoluteCenter /= N;
  
  target.resize( N );
  for ( size_t i=0; i < N; i++ ) { 
    target[i] = absolute_[i] - targetAbsoluteCenter;
  }

  isUpdatedE = false;
  isUpdatedGrad = false;
  isUpdated2ndDeriv = false;
  numberOfSavedRotations = 0;
}

void Fit::setRotor( const Vec3& rotor_ ) {
  if ( hasZeroRadius ) { return; }
  rotor = rotor_;
  double bmag2 = rotor.norm2();
  Vec3 tempRotor = rotor;
  
  //If the rotor's magnitude is greater than 2,
  //the rotor is invalid.  Here, we truncate
  //the rotor to magnitude 2, preserving its direction.
  //Notice that this is done on a temporary copy of the 
  //rotor.  The original rotor is left unchanged.
  if ( bmag2 > 4.0 ) {
    tempRotor *= 2.0/sqrt( bmag2 );
  }

  currentRotation.setRotor( tempRotor );
  size_t N = base.size();
  for ( size_t i = 0; i < N; i++ ) {
    currentRotation.rotate( base[i], rotated[i] );          
  }
  rotorIsZero = false;

  isUpdatedE = false;
  isUpdatedGrad = false;
  isUpdated2ndDeriv = false;
}

void Fit::collapseRotor() {
  if ( rotorIsZero ) return;

  if ( !numberOfSavedRotations++ ) {
    savedRotation.setRotor( currentRotation.rotor() );
  } 
  else {
    savedRotation.addRotation( currentRotation );
  }
  
  size_t N = base.size();
  for ( size_t i = 0; i < N; i++ ) {
    base[i] = rotated[i];          
  }
  rotor = Vec3( 0,0,0 );  
  rotorIsZero = true;
}

double Fit::energy() {
  if ( isUpdatedE ) return E;
  E = 0.0;
  size_t N = rotated.size();
  if ( target.size() != N ) {
    cout << "Error: Fit target size does not match source size" << endl;
    exit(0);
  }
  for ( size_t i = 0; i < N; i++ ) {
    E += rotated[i].dist2( target[i] );
  }
  E *= 0.5;
  
  //rotor stretch energy
  double b2 = rotor.norm2();
  double b;
  if ( b2 > 4.0 ) {
    b = sqrt(b2);
    E += 0.5*(b - 2.0)*(b - 2.0);
  }
  
  isUpdatedE = true;
  return E;
}

const Vec3& Fit::gradient() {
  if ( isUpdatedGrad ) return gradrotor;

  //handle special case - rigid unit has zero radius
  gradrotor = Vec3(0,0,0);
  if ( hasZeroRadius ) {
    isUpdatedGrad = true;
    return gradrotor;
  }

  //if rotor is big enough that it is in the unstable
  //region, collapse the rotor.  This also eliminates
  //the possibility that the rotor will have a magnitude
  //greater than 2 as we proceed with the gradient calculation
  if ( rotor.norm2() > rotormagcutoff2 ) {
    collapseRotor();
  }

  size_t N = rotated.size();
  if ( target.size() != N ) {
    cout << "Error: Fit target size does not match source size" << endl;
    exit(0);
  }
  Vec3 dV_dr;
  Rotator rotator;
  for ( size_t i = 0; i < N; i++ ) {
    dV_dr = rotated[i] - target[i];    
    // the gradient of the rotor x-component is
    // the dV_dx*dx_dbx + dV_dy*dy_dbx + dV_dz*dz_dbx.
    rotator.setRotor( rotor );
    gradrotor.x += 
      dV_dr.x * rotator.dx_dbx( base[i] ) +
      dV_dr.y * rotator.dy_dbx( base[i] ) +
      dV_dr.z * rotator.dz_dbx( base[i] );
    
    // the gradient of the rotor y and z components follow the
    // same pattern
    gradrotor.y += 
      dV_dr.x * rotator.dx_dby( base[i] ) +
      dV_dr.y * rotator.dy_dby( base[i] ) +
      dV_dr.z * rotator.dz_dby( base[i] );
    gradrotor.z += 
      dV_dr.x * rotator.dx_dbz( base[i] ) +
      dV_dr.y * rotator.dy_dbz( base[i] ) +
      dV_dr.z * rotator.dz_dbz( base[i] );
  }

  isUpdatedGrad = true;
  return gradrotor;
}

const Vec3& Fit::secondDerivativeDiagonal() {
  if ( isUpdated2ndDeriv ) return d2V_dB2;
  d2V_dB2 = Vec3(0,0,0);
  size_t N = rotated.size();
  Vec3 dV_dr;
  for ( size_t i = 0; i < N; i++ ) {
    dV_dr = rotated[i] - target[i];    
    d2V_dB2.x += 
      rotated[i].z*rotated[i].z +
      rotated[i].y*rotated[i].y -
      dV_dr.y*rotated[i].y - 
      dV_dr.z*rotated[i].z;
    d2V_dB2.y += 
      rotated[i].z*rotated[i].z +
      rotated[i].x*rotated[i].x -
      dV_dr.x*rotated[i].x - 
      dV_dr.z*rotated[i].z;
    d2V_dB2.z += 
      rotated[i].y*rotated[i].y +
      rotated[i].x*rotated[i].x -
      dV_dr.x*rotated[i].x - 
      dV_dr.y*rotated[i].y;
  }
  isUpdated2ndDeriv = true;
  return d2V_dB2;
}

void Fit::simpleFit() {
  savedRotation.setRotor( Vec3(0,0,0) );
  
  if ( hasZeroRadius ) {
    return;
  }
  
  collapseRotor();
  
  minimizer.enable_TolMaxPreconditionedGradComp( 0.001/radiusGyration );
  
  // a few initial steepest descent steps
  minimizer.setMethod_SteepestDescent();
  minimizer.minimize( func, 5 );
  int count = 0;
  // conjugate gradient steps
  minimizer.setMethod_ConjGrad();
  while ( !minimizer.isTolReached() && count < 400 ) {
    if ( count >= 50 && count%50==0 ) {
      minimizer.setMethod_SteepestDescent();
      minimizer.minimize( func, 5 );
      minimizer.setMethod_ConjGrad();
      collapseRotor();
    }
    minimizer.minimize( func, 3 );
    count++;
  }
  //cout << "Count " << count << endl;
  collapseRotor();

  if ( !minimizer.isTolReached() ) {
    Vec3 grad = gradient();
    cout << "Could not achieve desired gradient tolerance in fitting." << endl;
    cout << "Gradient " << grad << endl;
    for ( size_t i = 0; i < rotated.size(); i++ ) {
      cout << rotated[i] << " " << target[i] << endl;
    }
    exit(0);
  }
}

/*
void Fit::doFit() {
//do simple fit
//check 2nd derivative
//if negative, ping the rotor component that corresponds to
//  the negative term.
//repeat until simple fit returns with a positive definite
//  second derivative--WAIT some components may be zero,
//  and that should be ok (example- a line)
//
//
// Also - try testing the fit convergence with the preconditioner
// Is it better? 
}
*/
