#include "RigidUnitSystem.h"
#include "RigidUnitData.h"
#include "RigidUnitLookup.h"
#include <cmath>
#include <limits>
#include <cstdlib>

using namespace std;

RigidUnitSystem::RigidUnitSystem() :
  absolutePositionsChanged( false ),
  basePositionsChanged( false )
{
}

RigidUnitSystem::RigidUnitSystem(
            const vector< vector<int> >& RUtoPlist_,
            const vector<Vec3>& meanPositions_ ) :
  absolutePositionsChanged( false ),
  basePositionsChanged( false )
{
  reinitializeToNewMeanPositions( RUtoPlist_, meanPositions_ );
}            

RigidUnitSystem::~RigidUnitSystem()
{
}

void RigidUnitSystem::buildRigidUnitFromAbsolutePositions( int ru ) {

  const vector<int> *rup_list = 
    &lookup.getRUPlistFromRU( ru );

  int rup;
  if ( rup_list->size() == 1 ) {
    rup = (*rup_list)[0];
    data.hasZeroRadius[ru] = true;
    data.centers[ru] = 
      data.absolutePositions[rup];
    data.basePositions[rup] = Vec3(0.0,0.0,0.0);
    data.radius[ru] = 0;
    return;
  }
    
  Vec3 meanvec;
  meanvec = Vec3(0,0,0);
  for ( vector<int>::const_iterator rup_it = 
    rup_list->begin();
  rup_it != rup_list->end();
  rup_it++ ) {
    meanvec += data.absolutePositions[*rup_it];    
  }
  meanvec /= static_cast<double>( rup_list->size() );
  data.centers[ru] = meanvec;

  double max_r2 = 0.0;
  double r2;
  for ( vector<int>::const_iterator rup_it = 
    rup_list->begin();
  rup_it != rup_list->end();
  rup_it++ ) {
    rup = *rup_it;
    data.basePositions[rup] = 
      data.absolutePositions[rup] - meanvec;
    r2 = data.basePositions[rup].norm2();
    if ( r2 > max_r2 ) max_r2 = r2;
  }
  data.radius[ru] = sqrt(max_r2);
  data.hasZeroRadius[ru] = 
    data.radius[ru] <= numeric_limits<double>::epsilon();
}

void RigidUnitSystem::reinitializeToNewMeanPositions(
           const vector< vector<int> >& RUtoPlist_,
           const vector<Vec3>& meanPositions_ )
{  
  size_t nP = meanPositions_.size();
  lookup.setLookup_RUtoPlist( RUtoPlist_, nP );
  int nRUP = lookup.nRigidUnitPoints();
  int nRU = lookup.nRigidUnits();
  
  //initialize mean positions
  data.meanPositions = meanPositions_;
  
  //initialize rigid unit point absolute positions
  data.absolutePositions.resize( nRUP );
  for ( int rup = 0; rup < nRUP; rup++ ) {
    int p = lookup.getPfromRUP( rup );
    data.absolutePositions[rup] = data.meanPositions[p];
  }


  data.rotors.resize( nRU, Vec3(0,0,0) );
  data.centers.resize( nRU );
  data.basePositions.resize( nRUP );
  data.radius.resize( nRU );
  data.hasZeroRadius.resize( nRU );
 
  for ( int ru = 0; ru < nRU; ru++ ) {
    buildRigidUnitFromAbsolutePositions( ru );
  }

  absolutePositionsChanged = true;
  basePositionsChanged = true;
}

void RigidUnitSystem::collapseRotor( int ru ) {
  if ( abs(data.rotors[ru].x) < numeric_limits<double>::epsilon() &&
       abs(data.rotors[ru].y) < numeric_limits<double>::epsilon() &&
       abs(data.rotors[ru].z) < numeric_limits<double>::epsilon() ) return;
       
  const vector<int> *rup_list =
    &lookup.getRUPlistFromRU( ru );
  int rup;
  for ( vector<int>::const_iterator rup_it = 
          rup_list->begin(); rup_it != rup_list->end(); rup_it++ ) {
    rup = *rup_it;
    data.basePositions[rup] = 
      data.absolutePositions[rup] - data.centers[ru];
  }
  data.rotors[ru] = Vec3(0,0,0);
  basePositionsChanged = true;
}

void RigidUnitSystem::collapseRotors() {
  int nRU = lookup.nRigidUnits();
  for ( int ru = 0; ru < nRU; ru++ ) {
    collapseRotor( ru );
  }
}

void RigidUnitSystem::updateAbsolutePositions( int ru ) {
  Vec3 tempRotor = data.rotors[ru];
  double bmag2 = tempRotor.norm2();
  
  //If the rotor's magnitude is greater than 2,
  //the rotor is invalid.  Here, we truncate
  //the rotor to magnitude 2, preserving its direction.
  //Notice that this is done on a temporary copy of the 
  //rotor.  The original rotor is left unchanged.
  if ( bmag2 > 4.0 ) {
    tempRotor *= 2.0/sqrt( bmag2 );
  }

  Rotator rotator( tempRotor );

  const vector<int> *rup_list =
    &lookup.getRUPlistFromRU( ru );
  int rup;
  for ( vector<int>::const_iterator rup_it = 
    rup_list->begin(); rup_it != rup_list->end(); rup_it++ ) {
    rup = *rup_it; 
    rotator.rotate( data.basePositions[rup],
                    data.absolutePositions[rup] );          
    data.absolutePositions[rup] += data.centers[ru];
  }
  absolutePositionsChanged = true;
}

void RigidUnitSystem::updateAbsolutePositions() {
  int nRU = lookup.nRigidUnits();
  #pragma omp parallel for
  for ( int ru = 0; ru < nRU; ru++ ) {
    updateAbsolutePositions( ru );
  }
}

void RigidUnitSystem::setCenters( 
               const vector<Vec3> &centers_ ) {
  if ( centers_.size() != lookup.nRigidUnits() ) {
    cout << "Error: Number of input vectors \n";
    cout << "  does not match the number of rigid units in \n";
    cout << "  the AbstractRigidUnitLookup object" << endl;
    exit(0);
  }
  data.centers = centers_;
  updateAbsolutePositions();
}

void RigidUnitSystem::setRotors( 
               const vector<Vec3> &rotors_ ) {
  if ( rotors_.size() != lookup.nRigidUnits() ) {
    cout << "Error: Number of input vectors \n";
    cout << "  does not match the number of rigid units in \n";
    cout << "  the AbstractRigidUnitLookup object" << endl;
    exit(0);
  }
  data.rotors = rotors_;
  updateAbsolutePositions();
}

void RigidUnitSystem::setCentersAndRotors( 
      const vector<Vec3> &centers_, const vector<Vec3> &rotors_ ) {
  if ( centers_.size() != lookup.nRigidUnits() ||
       rotors_.size() != lookup.nRigidUnits() ) {  
    cout << "Error: Number of input vectors \n";
    cout << "  does not match the number of rigid units in \n";
    cout << "  the AbstractRigidUnitLookup object" << endl;
    exit(0);
  }
  data.centers = centers_;
  data.rotors = rotors_;
  updateAbsolutePositions();
}

void RigidUnitSystem::setCenter( int ru, const Vec3& center_ ) {
  data.centers[ru] = center_;
  updateAbsolutePositions( ru );
}
void RigidUnitSystem::setRotor( int ru, const Vec3& rotor_ ) {
  data.rotors[ru] = rotor_;
  updateAbsolutePositions( ru );
}
void RigidUnitSystem::setCenterAndRotor( int ru,
             const Vec3& center_,
             const Vec3& rotor_ ) {
  data.centers[ru] = center_;
  data.rotors[ru] = rotor_;
  updateAbsolutePositions( ru );
}

void RigidUnitSystem::addToRotor( int ru, 
            const Vec3& deltaRotor ) {
  data.rotors[ru] += deltaRotor;
  updateAbsolutePositions( ru );
}

void RigidUnitSystem::addToCenter( int ru,
             const Vec3& deltaCenter ) {
  data.centers[ru] += deltaCenter;
  updateAbsolutePositions( ru );
}


void RigidUnitSystem::clearState() {
  absolutePositionsChanged = false;
  basePositionsChanged = false;
}

void RigidUnitSystem::rotate( int ru,
              Rotator& rotator,
              const Vec3& centerOfRotation ) {

  //If rotor is not zero, collapse the rotor.
  if ( abs( data.rotors[ru].x ) > numeric_limits<double>::epsilon() ||
       abs( data.rotors[ru].y ) > numeric_limits<double>::epsilon() ||
       abs( data.rotors[ru].z ) > numeric_limits<double>::epsilon() ) {
    collapseRotor( ru );
  }

  //now calculate and apply the "orbital" rotation:
  //the rigid unit's center rotated about the input center of rotation.
  Vec3 relativeCenter = data.centers[ru] - centerOfRotation;
  Vec3 newRelativeCenter;
  rotator.rotate( relativeCenter, newRelativeCenter );
  data.centers[ru] = newRelativeCenter + centerOfRotation;

  //now apply the "spin" rotation: the rigid unit
  //about its center.
  data.rotors[ru] = rotator.rotor();
  const vector<int> *rup_list =
    &lookup.getRUPlistFromRU( ru );
  int rup;
  for ( vector<int>::const_iterator rup_it = 
          rup_list->begin(); rup_it != rup_list->end(); rup_it++ ) {
    rup = *rup_it;
    rotator.rotate( data.basePositions[rup],
                    data.absolutePositions[rup] );          
    data.absolutePositions[rup] += data.centers[ru];
  }
  absolutePositionsChanged = true;

  //To finish, collapse the rotor to zero.
  collapseRotor( ru );

}

void RigidUnitSystem::updateMeanPositions()
{
  
  size_t nP = lookup.nPoints();
  #pragma omp parallel for
  for ( size_t p = 0; p < nP; p++ ) {
    const vector<int> *rup_list = 
      &lookup.getRUPlistFromP( p );
    vector<int>::const_iterator rup_it =
       rup_list->begin();
    int rup = *rup_it;
    data.meanPositions[p] = data.absolutePositions[rup];

    if ( ++rup_it != rup_list->end() ) {
      //we already have accounted for the first element
      //so here we begin with next element.
      for ( ; rup_it != rup_list->end(); rup_it++ ) {
        rup = *rup_it;
        data.meanPositions[p] += data.absolutePositions[rup];
      } 
      data.meanPositions[p] /= static_cast<double>( rup_list->size() );
    }
  }  
}

