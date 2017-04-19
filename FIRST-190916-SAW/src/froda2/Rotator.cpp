#include "Rotator.h"
#include <cmath>
#include <limits>
#include <cstdlib>

using namespace std;

Rotator::Rotator() : 
  bmag2_lim( 4.0*(1.0 + 100.0*numeric_limits<double>::epsilon() ) ),
  isRotationMatrixUpdated( false ),
  areUpdatedFirstDerivatives( false )  
{}

Rotator::Rotator( const Vec3& rotor_ ) :
  bmag2_lim( 4.0*(1.0 + 10.0*numeric_limits<double>::epsilon() ) ),
  isRotationMatrixUpdated( false ),
  areUpdatedFirstDerivatives( false )  
{
  setRotor( rotor_ );
}

Rotator::~Rotator()
{
}

void Rotator::setRotor( const Vec3& rotor_ ) {
  b = rotor_;
  bxbx = b.x*b.x;
  byby = b.y*b.y;
  bzbz = b.z*b.z;
  bxby = b.x*b.y;
  bybz = b.y*b.z;
  bzbx = b.z*b.x;
  double bmag2 = bxbx + byby + bzbz;
  
  if ( bmag2 >= bmag2_lim ) {
    cerr << "Rotor magnitude is greater than 2.0. " << endl;
    cerr << sqrt( bmag2 ) << endl;
    exit(0);
  }
  
  double X2 = 1.0 - 0.25 * bmag2;

  // if we reach here, we know that the rotor magnitude
  // is <= 2.0.  However, if the magnitude is within machine
  // precision of 2.0, it is possible that the quantity X-squared
  // is calculated to be slightly negative instead of zero.
  // Taking the sqrt of a negative number will generate NaN,
  // so to prevent that we just set X to zero.
  if ( X2 <= 0.0 ) {
    X = 0;
  }
  else {
    X = sqrt( X2 );
  }

  if ( isnan(X) ) {
    cerr << "X is nan" << endl;
    cerr << "Rotor is " << b << endl;
    exit(0);
  }
  
  updateRotationMatrix();
  
  areUpdatedFirstDerivatives = false;
}

void Rotator::updateRotationMatrix() {
  Vec3 b_TimesX = b * X;
  
  M11 = 1.0 - 0.5*(byby + bzbz);
  M12 = -b_TimesX.z + 0.5*bxby;
  M13 = b_TimesX.y + 0.5*bzbx;

  M21 = b_TimesX.z + 0.5*bxby; 
  M22 = 1.0 - 0.5*(bxbx + bzbz); 
  M23 = -b_TimesX.x + 0.5*bybz;

  M31 = -b_TimesX.y + 0.5*bzbx; 
  M32 = b_TimesX.x + 0.5*bybz; 
  M33 = 1.0 - 0.5*(bxbx + byby);
  
  isRotationMatrixUpdated = true;
}

void Rotator::updateFirstDerivatives() {
  if ( X <= numeric_limits<double>::epsilon() ) {
    cout << "Error: Cannot take derivative with repect to rotor parameters" << endl;
    cout << "       when rotor mag >= 2.0." << endl;
    exit(0);
  }
  double X_times_4 = X*4.0;
  double bx_over_2 = b.x/2.0;
  double by_over_2 = b.y/2.0;
  double bz_over_2 = b.z/2.0;
  double bxbx_over_4X = bxbx/X_times_4;
  double byby_over_4X = byby/X_times_4;
  double bzbz_over_4X = bzbz/X_times_4;
  double bxby_over_4X = bxby/X_times_4;
  double bybz_over_4X = bybz/X_times_4;
  double bzbx_over_4X = bzbx/X_times_4;
  
  dx_dbx_multiplier.x = 0.0;
  dx_dbx_multiplier.y = bzbx_over_4X + by_over_2;
  dx_dbx_multiplier.z = -bxby_over_4X + bz_over_2;
  
  dx_dby_multiplier.x = -b.y;
  dx_dby_multiplier.y = bybz_over_4X + bx_over_2;
  dx_dby_multiplier.z = X - byby_over_4X;
  
  dx_dbz_multiplier.x = -b.z;
  dx_dbz_multiplier.y = -X + bzbz_over_4X;
  dx_dbz_multiplier.z = -bybz_over_4X + bx_over_2;
  
  dy_dbx_multiplier.x = -bzbx_over_4X + by_over_2;
  dy_dbx_multiplier.y = -b.x;
  dy_dbx_multiplier.z = -X + bxbx_over_4X;
  
  dy_dby_multiplier.x = -bybz_over_4X + bx_over_2;
  dy_dby_multiplier.y = 0.0;
  dy_dby_multiplier.z = bxby_over_4X + bz_over_2;
  
  dy_dbz_multiplier.x = X - bzbz_over_4X;
  dy_dbz_multiplier.y = -b.z;
  dy_dbz_multiplier.z = bzbx_over_4X + by_over_2;
  
  dz_dbx_multiplier.x = bxby_over_4X + bz_over_2;
  dz_dbx_multiplier.y = X - bxbx_over_4X;
  dz_dbx_multiplier.z = -b.x;
  
  dz_dby_multiplier.x = -X + byby_over_4X;
  dz_dby_multiplier.y = -bxby_over_4X + bz_over_2;
  dz_dby_multiplier.z = -b.y;
  
  dz_dbz_multiplier.x = bybz_over_4X + bx_over_2;
  dz_dbz_multiplier.y = -bzbx_over_4X + by_over_2;
  dz_dbz_multiplier.z = 0.0;
  
  areUpdatedFirstDerivatives = true;
} 

void Rotator::addRotation( Rotator& other ) {
  
  if ( !isRotationMatrixUpdated ) {
    updateRotationMatrix();
  }
  if ( !other.isRotationMatrixUpdated ) {
    other.updateRotationMatrix();
  }
  
  //multiply matrices
  double C11 = other.M11*M11 + other.M12*M21 + other.M13*M31;
  double C12 = other.M11*M12 + other.M12*M22 + other.M13*M32;
  double C13 = other.M11*M13 + other.M12*M23 + other.M13*M33;
  double C21 = other.M21*M11 + other.M22*M21 + other.M23*M31;
  double C22 = other.M21*M12 + other.M22*M22 + other.M23*M32;
  double C23 = other.M21*M13 + other.M22*M23 + other.M23*M33;
  double C31 = other.M31*M11 + other.M32*M21 + other.M33*M31;
  double C32 = other.M31*M12 + other.M32*M22 + other.M33*M32;
  double C33 = other.M31*M13 + other.M32*M23 + other.M33*M33;
    
  //construct new rotor
  double X_times_2_squared = 1.0 + C11 + C22 + C33;
  Vec3 rotorEquiv;
  if ( X_times_2_squared > 100*numeric_limits<double>::epsilon() ) {
    double X_times_2 = sqrt( X_times_2_squared );
    rotorEquiv.x = ( C32 - C23 )/X_times_2;
    rotorEquiv.y = ( C13 - C31 )/X_times_2;
    rotorEquiv.z = ( C21 - C12 )/X_times_2;
  }
  else {
    //this is a special case when the equivalent rotation is very close to 
    //180 degrees.  In this case, X is near zero, so division by X 
    //as performed above would yield
    //infinity.  Also, in this case, we know that the square magnitude of 
    //b (the rotor) is about 4.
    //
    //First, we calculate the square of the rotor components, bx2, by2, 
    //and bz2, in a manner that does not
    //risk division by zero, looking for any one of them that is larger
    //than some threshold.  We will use this component as a reference
    //to determine the other two components.  For example, say that
    //we find that bx squared is greater than some threshold.  Then
    //we will use bx to determine by and bz.  We arbitrarily choose
    //the positive root of bx squared as the value of bx.
    
    //a note on the threshold--because we know that b2 = 4,
    //one of the square components of b must be >= 4/3.
    //Here, we set the threshold to 1, as there must certainly
    //be a component whose square magnitude is greater than 1.
    double threshold = 1.0;
    
    double bx2 = ( C11 + 1 )*2.0;
    if ( bx2 > threshold ) {
      rotorEquiv.x = sqrt( bx2 );
      rotorEquiv.y = ( C12 + C21 )/rotorEquiv.x;
      rotorEquiv.z = ( C13 + C31 )/rotorEquiv.x;
    }
    else {
      double by2 = ( C22 + 1 )*2.0;
      if ( by2 > threshold ) {
        rotorEquiv.y = sqrt( by2 );
        rotorEquiv.x = ( C12 + C21 )/rotorEquiv.y;
        rotorEquiv.z = ( C23 + C32 )/rotorEquiv.y;
      }
      else {
        double bz2 = ( C33 + 1 )*2.0;
        if ( bz2 > threshold ) {
          rotorEquiv.z = sqrt( bz2 );
          rotorEquiv.x = ( C13 + C31 )/rotorEquiv.z;
          rotorEquiv.y = ( C32 + C23 )/rotorEquiv.z;
        }
        else {
          cout << "Rotator addition error" << endl;
          exit(0);
        }
      }
    }
  
  }

  double mag2 = rotorEquiv.norm2();
  if ( mag2 > 4.0 ) {
    rotorEquiv *= 2.0/sqrt( mag2 );
  }
    
  //update internal rotor
  setRotor( rotorEquiv );

  //update internal rotation matrix
  M11 = C11;
  M12 = C12;
  M13 = C13;
  M21 = C21;
  M22 = C22;
  M23 = C23;
  M31 = C31;
  M32 = C32;
  M33 = C33;  
  isRotationMatrixUpdated = true;
  
  //update other member variables
  areUpdatedFirstDerivatives = false;
}

