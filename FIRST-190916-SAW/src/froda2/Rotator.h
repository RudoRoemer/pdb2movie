#ifndef ROTATOR_H_
#define ROTATOR_H_

#include "Vec3.h"

class Rotator;

class Rotator
{
public:
  Rotator();
	Rotator( const Vec3& rotor_ );
	virtual ~Rotator();
  
  void setRotor( const Vec3& rotor_ );
  const Vec3& rotor() const { return b; } 

  void addRotation( Rotator& other );

  void rotate( const Vec3& baseVec, Vec3& rotatedVec ) {
    if ( !isRotationMatrixUpdated ) {
      updateRotationMatrix();
    }
    rotatedVec.x = M11*baseVec.x + M12*baseVec.y + M13*baseVec.z;
    rotatedVec.y = M21*baseVec.x + M22*baseVec.y + M23*baseVec.z;
    rotatedVec.z = M31*baseVec.x + M32*baseVec.y + M33*baseVec.z;
  }

  double dx_dbx( const Vec3& base ) {
    if ( !areUpdatedFirstDerivatives ) {
      updateFirstDerivatives();
    }
    return dx_dbx_multiplier.dot( base );
  }
  double dx_dby( const Vec3& base ) {
    if ( !areUpdatedFirstDerivatives ) {
      updateFirstDerivatives();
    }
    return dx_dby_multiplier.dot( base );
  }
  double dx_dbz( const Vec3& base ) {
    if ( !areUpdatedFirstDerivatives ) {
      updateFirstDerivatives();
    }
    return dx_dbz_multiplier.dot( base );
  }
  double dy_dbx( const Vec3& base ) {
    if ( !areUpdatedFirstDerivatives ) {
      updateFirstDerivatives();
    }
    return dy_dbx_multiplier.dot( base );
  }
  double dy_dby( const Vec3& base ) {
    if ( !areUpdatedFirstDerivatives ) {
      updateFirstDerivatives();
    }
    return dy_dby_multiplier.dot( base );
  }
  double dy_dbz( const Vec3& base ) {
    if ( !areUpdatedFirstDerivatives ) {
      updateFirstDerivatives();
    }
    return dy_dbz_multiplier.dot( base );
  }
  double dz_dbx( const Vec3& base ) {
    if ( !areUpdatedFirstDerivatives ) {
      updateFirstDerivatives();
    }
    return dz_dbx_multiplier.dot( base );
  }
  double dz_dby( const Vec3& base ) {
    if ( !areUpdatedFirstDerivatives ) {
      updateFirstDerivatives();
    }
    return dz_dby_multiplier.dot( base );
  }
  double dz_dbz( const Vec3& base ) {
    if ( !areUpdatedFirstDerivatives ) {
      updateFirstDerivatives();
    }
    return dz_dbz_multiplier.dot( base );
  }

private:
  Vec3 b;
  double X;
  double bxbx;
  double byby;
  double bzbz;
  double bxby;
  double bybz;
  double bzbx;
  double M11;
  double M12;
  double M13;
  double M21;
  double M22;
  double M23;
  double M31;
  double M32;
  double M33;
  Vec3 dx_dbx_multiplier;
  Vec3 dx_dby_multiplier;
  Vec3 dx_dbz_multiplier;
  Vec3 dy_dbx_multiplier;
  Vec3 dy_dby_multiplier;
  Vec3 dy_dbz_multiplier;
  Vec3 dz_dbx_multiplier;
  Vec3 dz_dby_multiplier;
  Vec3 dz_dbz_multiplier;
  const double bmag2_lim;
  bool isRotationMatrixUpdated;
  bool areUpdatedFirstDerivatives;
  void updateRotationMatrix();
  void updateFirstDerivatives();
};

#endif /*ROTATOR_H_*/
