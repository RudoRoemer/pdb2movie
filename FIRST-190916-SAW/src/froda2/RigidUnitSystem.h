#ifndef RIGIDUNITSYSTEM_H_
#define RIGIDUNITSYSTEM_H_

#include "RigidUnitData.h"
#include "RigidUnitLookup.h"
#include "Observable.h"
#include "Rotator.h"
#include <vector>
#include "Vec3.h"

class RigidUnitSystem : public Observable
{
public:
  RigidUnitSystem(); 

	RigidUnitSystem( 
           const std::vector< std::vector<int> >& RUtoPlist_,
           const std::vector<Vec3>& meanPositions_ );
           
	virtual ~RigidUnitSystem();
  
  const std::vector<Vec3>& meanPositions() const {
    return data.meanPositions;
  }
  const Vec3& meanPositions( int p ) const {
    return data.meanPositions[p];
  }
  const std::vector<Vec3>& rotors() const {
    return data.rotors;
  }
  const Vec3& rotors( int ru ) const {
    return data.rotors[ru];
  }
  const std::vector<Vec3>& centers() const {
    return data.centers;
  }
  const Vec3& centers( int ru ) const {
    return data.centers[ru];
  }
  const std::vector<Vec3>& basePositions() const {
    return data.basePositions;
  }
  const Vec3& basePositions( int rup ) const {
    return data.basePositions[rup];
  }
  const std::vector<Vec3>& absolutePositions() const {
    return data.absolutePositions;
  }
  const Vec3& absolutePositions( int rup ) const {
    return data.absolutePositions[rup];
  }
  double radius( int ru ) const {
    return data.radius[ru];
  }
  bool hasZeroRadius( int ru ) const {
    return static_cast<bool>( data.hasZeroRadius[ru] );
  }
  
  const std::vector< int >& getPlistFromRU( int ru ) const {
    return lookup.getPlistFromRU( ru );
  }
  const std::vector< int >& getRUlistFromP( int p ) const {
    return lookup.getRUlistFromP( p );
  }
  const std::vector< int >& getRUPlistFromP( int p ) const {
    return lookup.getRUPlistFromP( p );
  }
  const int getPfromRUP( int rup ) const {
    return lookup.getPfromRUP( rup );
  }
  const size_t nPoints() const {
    return lookup.nPoints();
  }
  const std::vector< int >& getRUPlistFromRU( int ru ) const {
    return lookup.getRUPlistFromRU( ru );
  }
  int getRUfromRUP( int rup ) const {
    return lookup.getRUfromRUP( rup );
  }
  size_t nRigidUnitPoints() const {
    return lookup.nRigidUnitPoints();
  }
  size_t nRigidUnits() const {
    return lookup.nRigidUnits();
  }  
  bool doPointsBelongToSameRigidUnit( int p1, int p2 ) const {
    return lookup.doPointsBelongToSameRigidUnit( p1, p2 );
  }


  bool AbsolutePositionsChanged() const { 
    return absolutePositionsChanged;
  }
  bool BasePositionsChanged() const { 
    return basePositionsChanged;
  }

  void collapseRotors();
    
  void setCenters( const std::vector<Vec3> &centers_ );
  void setRotors( const std::vector<Vec3> &rotors_ );
  void setCentersAndRotors( const std::vector<Vec3> &centers_, 
                            const std::vector<Vec3> &rotors_ );
  void setCenter( int ru, const Vec3& center_ );
  void setRotor( int ru, const Vec3& rotor_ );
  void setCenterAndRotor( int ru,
        const Vec3& center_,
        const Vec3& rotor_ );
  void addToRotor( int ru, const Vec3& deltaRotor );
  void addToCenter( int ru, const Vec3& deltaCenter );  
  void rotate( int ru,
         Rotator& rotator,
         const Vec3& absoluteCenterOfRotation );
  void collapseRotor( int ru );
  void clearState();

  void update() {
    updateMeanPositions();
    notifyObservers();
    clearState();
  }
  
  void reinitializeToNewMeanPositions( 
           const std::vector< std::vector<int> >& RUtoPlist_,
           const std::vector<Vec3>& meanPositions_ );  


private:
  RigidUnitData data;
  RigidUnitLookup lookup;

  bool absolutePositionsChanged;
  bool basePositionsChanged;
  
  void updateAbsolutePositions( int ru );
  void updateAbsolutePositions();
  void buildRigidUnitFromAbsolutePositions( int ru ); 
  void updateMeanPositions();
};

#endif /*RIGIDUNITSYSTEM_H_*/
