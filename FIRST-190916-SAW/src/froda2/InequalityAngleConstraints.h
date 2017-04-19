#ifndef INEQUALITYANGLECONSTRAINTS_H_
#define INEQUALITYANGLECONSTRAINTS_H_

#include <vector>
#include "Energy.h"
#include "Gradient.h"
#include "Mismatch.h"
#include "Vec3.h"
#include "GeneralizedCoords.h"
#include <cmath>
class RigidUnitSystem;

class InequalityAngleConstraint {
public:
  InequalityAngleConstraint( int p1_, int p2_, int p3_, double minAllowedTheta_ ) :
    p1(p1_), p2(p2_), p3(p3_), 
    minAllowedTheta(minAllowedTheta_), 
    maxAllowedCosTheta(cos(minAllowedTheta_)) {}
  InequalityAngleConstraint& operator=( const InequalityAngleConstraint& other ) {
    p1 = other.p1;
    p2 = other.p2;
    p3 = other.p3;
    minAllowedTheta = other.minAllowedTheta;
    maxAllowedCosTheta = other.maxAllowedCosTheta;
    return *this;
  }
  int p1;
  int p2;
  int p3;
  double minAllowedTheta;
  double maxAllowedCosTheta;
};

class InequalityAngleConstraints : 
  public EnergyTerm, 
  public GradientTerm,
  public MismatchTerm
{
public:
	InequalityAngleConstraints( RigidUnitSystem *rigidUnitSystem );
	virtual ~InequalityAngleConstraints();

	double energy();
  void addToGradient( std::vector<Vec3> &dV_dr_rigidUnitPoint,
                      std::vector<SecondDerivative> &secondDerivative_rigidUnitPoint );
  double mismatch();

  void insert( const InequalityAngleConstraint& con ) {
    angles.push_back( con );
  }

  void insert( int p1, int p2, int p3, double minAllowedTheta ) {
    angles.push_back( InequalityAngleConstraint( p1, p2, p3, minAllowedTheta ) );
  }
  
  void insert( const InequalityAngleConstraints& otherInequalityAngleConstraints ) {
    angles.insert( angles.end(),
        otherInequalityAngleConstraints.begin(),
        otherInequalityAngleConstraints.end() );
  }
  
  typedef std::vector<InequalityAngleConstraint>::iterator iterator;
  typedef std::vector<InequalityAngleConstraint>::const_iterator const_iterator;

  iterator begin() { return angles.begin(); }
  const_iterator begin() const { return angles.begin(); }
  iterator end() { return angles.end(); }
  const_iterator end() const { return angles.end(); }
  size_t size() const { return angles.size(); }
  InequalityAngleConstraint& operator[]( size_t i ) { return angles[i]; }
  const InequalityAngleConstraint& operator[]( size_t i ) const { return angles[i]; }

private:
  std::vector<InequalityAngleConstraint> angles;
  const RigidUnitSystem *rigidUnitSystem;

};

#endif /*INEQUALITYANGLECONSTRAINTS_H_*/
