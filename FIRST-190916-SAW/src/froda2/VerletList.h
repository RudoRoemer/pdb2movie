#ifndef _VERLETLIST_H_
#define _VERLETLIST_H_

#include "MovingPoint.h"
#include "Spanner.h"
#include "Vec3.h"
#include <cmath>

// class VerletListPoint derives from class MovingPoint.
// It merely adds an interaction cutoff distance parameter,
// used in the VerletList class.
//////////////////////////////////////////////////////////////////////
class VerletListPoint : public MovingPoint {
public:
  VerletListPoint() : MovingPoint(), cutoffDist(0) {}
  VerletListPoint( const Vec3 *position, double cutoffDist, int id );
  VerletListPoint( const VerletListPoint& other );
  ~VerletListPoint();
  const VerletListPoint& operator=( const VerletListPoint& other );
  double getCutoffDist() const {return cutoffDist;}
private:
  double cutoffDist;
};

class VL_NeighborInfo {
 public:
  int id;
  double pairCutoff;
  //double deltaSquared;
  //Vec3 delta21;
};

// class VerletList is used for keeping track of
// pairs of points that are within interacting distance.  It uses
// a smart updating feature that only refreshes its internal lists of
// interacting points when necessary, described below.  The user interacts
// with the VerletList as follows:
//   First, insert points into the VerletList, specifying the interaction
//     cutoff distance of each.
//   Commit the VerletList, which locks the list of points and puts the
//     points into a spanner.  The user here must supply the desired
//     spanner parameters, and also the desired cushion parameter (see below).
//   Update the proximity monitor as points move.
//   Query the VerletList for the set of points that are within the
//     interaction cutoff distance of a particular point.  You will get some
//     extras here, particularly points that were in the "cushion region"
//     (see below).
//
// Here's how the smart updating works.  These details are abstracted away from the user.
// When the user calls the update
// member function, the VerletList queries the Spanner for a list of
// nearby points relative to every point.  But it does not just query out to the
// interaction cutoff of each point.  Instead, it queries out to the cutoff plus
// an extra cushion distance.  For each point, the VerletList stores a list of nearby points,
// which includes any extra points that lie in the cushion region.  As long as points
// do not move more than cushion/2, this stored list of nearby points is guaranteed
// to contain all points within the cutoff distance.
//
// When the user subsequently calls update(), the VerletList will check 
// to see how far points have moved.  If points have not moved more than cushion/2,
// then it will not even bother querying the spanner.  These update() calls can be
// very fast.  As points keep moving, at some point one of the points will wander
// more than cutoff/2.  Then, when update() is called, the VerletList knows that
// its lists of nearby points are no longer guaranteed to be accuate, so it re-queries
// the spanner to get up-to-date lists.
//
////////////////////////////////////////////////////////////////////////////////
template < class _IsPairExcluded, class _PairInteractionCutoff >
class VerletList {

public:
  typedef vector< vector<VL_NeighborInfo> > ProximityLists;
  
  VerletList( _IsPairExcluded &isPairExcluded_,
              _PairInteractionCutoff &pairInteractionCutoff_ );
  virtual ~VerletList();

  const vector<VL_NeighborInfo>& getProximityList( int point1 ) const;
  void update();
  void free();
  void commit( double beta_, double c_, 
               double interactionCutoff_, double cushion_ );
  void insert( const Vec3 *position, double cutoffDist, int id );
  void verify() const;
  void printstat() const;

private:
  Spanner *spanner;
  map<int,int> mapExternalIDtoInternalID;
  ProximityLists proximityLists;
  vector<VerletListPoint> verletListPoints;
  _IsPairExcluded &isPairExcluded;
  _PairInteractionCutoff &pairInteractionCutoff;
  
  double cushion;
  float beta;
  float c;
  float unitDistance;
  float movementThreshhold_VerletList;
  float movementThreshhold_Spanner;
  const vector<VL_NeighborInfo> emptylist;
  
  void updatePastPositions();
  double calcMaxMoved() const;
  void updateProximityLists();
};

////////////////////////////////////////////////////////////////////////////////
// Description:
//   VerletListPoint Constructor.  
// Parameters:
//   position - the address of the position 3-D Vector
//   cutoffDist - the center-to-center distance from the point beyond which
//     no interactions occur.  
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
inline VerletListPoint::VerletListPoint( const Vec3 *position_, double cutoffDist_, int id_ ) : 
  MovingPoint( position_, id_ ),
  cutoffDist(cutoffDist_) {}

inline VerletListPoint::VerletListPoint( const VerletListPoint& other ) : MovingPoint(other),
     cutoffDist(other.cutoffDist) {}
////////////////////////////////////////////////////////////////////////////////
// Description:
//   defines the '=' operator for VerletListPoints.
////////////////////////////////////////////////////////////////////////////////
inline const VerletListPoint& VerletListPoint::operator=(const VerletListPoint& other) {
  position = other.position;
  pastPosition = other.pastPosition;
  id = other.id;
  cutoffDist = other.cutoffDist;
  return other;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   VerletListPoint Destructor.  
////////////////////////////////////////////////////////////////////////////////
inline VerletListPoint::~VerletListPoint() {}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Get the list of points that are within the cutoff distance of a particular point
// Parameters:
//   point1 - id of particular point
// Return Value List:
//   returns all points that are within the cutoff distance of the input point.
//   Can also include extra points that are in the cushion region, outside the
//   cutoff distance.
////////////////////////////////////////////////////////////////////////////////
template < class A, class B >
inline const vector<VL_NeighborInfo>& VerletList<A,B>::getProximityList( int point1 ) const {
  map<int,int>::const_iterator iter = mapExternalIDtoInternalID.find( point1 );
  if ( iter == mapExternalIDtoInternalID.end() ) {
    return emptylist;
  }
  return proximityLists[iter->second];
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   VerletList Constructor.  
// Parameters:
//   none
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
template < class _IsPairExcluded, class _PairInteractionCutoff >
inline VerletList<_IsPairExcluded, _PairInteractionCutoff>::VerletList( 
    _IsPairExcluded &isPairExcluded_,
    _PairInteractionCutoff &pairInteractionCutoff_ ) : 
                spanner(NULL),
                isPairExcluded(isPairExcluded_),
                pairInteractionCutoff(pairInteractionCutoff_),
                cushion(0),
					      beta(0),
					      c(0),
					      unitDistance(0),
                movementThreshhold_VerletList(0.0),
                movementThreshhold_Spanner(0.0),
					      emptylist(0) {}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   VerletListPoint Destructor.  
////////////////////////////////////////////////////////////////////////////////
template < class A, class B >
inline VerletList<A, B>::~VerletList<A, B>() {
  delete spanner;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Insert a point into the VerletList.  Points can be inserted only when
//   setting up the VerletList.  After you call commit(), you cannot
//   insert points anymore.
// Parameters:
//   position - the address of the point's 3-D Vector.  As the point moves,
//     the proximity monitor will always know what its new position is.
//     IMPORTANT - this address must always be valid!  If the 3-D Vector
//     ends up moving to a different location in memory, this will be bad.
//     For example, if the 3-D Vector is stored in a std::vector<>, and then
//     you extend the size of the std::vector, the std::vector may move its
//     contents somewhere else, invalidating all pointers to its contents.
//   cutoffDist - defines the center-to-center distance beyond which points are
//     not neighbors.  When you query the VerletList for this point, it will return points
//     that are within cutoffDist of this point.
//   id - the ID of the point
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
template < class A, class B >
inline void VerletList<A, B>::insert( const Vec3 *position, double cutoffDist, int externalID ) {

  if ( spanner != NULL ) {
    cerr << "Error, cannot insert point into committed proximity monitor" << endl;
    exit(0);
  }
  
  map<int,int>::const_iterator iter = mapExternalIDtoInternalID.find( externalID );
  if ( iter != mapExternalIDtoInternalID.end() ) {
    cerr << "Error, attempting to insert duplicate point into proximity monitor " << endl;
    exit(0);
  }
  int internalID = verletListPoints.size();
  mapExternalIDtoInternalID[externalID] = internalID;
  verletListPoints.push_back( VerletListPoint( position, cutoffDist, externalID ) );

}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Clears the internal spanner, the stored points, and the stored
//   proximity lists
// Parameters:
//   none
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
template < class _IsPairExcluded, class _PairInteractionCutoff >
inline void VerletList<_IsPairExcluded, _PairInteractionCutoff>::free() {
  delete spanner;
  spanner = NULL;
  mapExternalIDtoInternalID.clear();
  verletListPoints.clear();
  proximityLists.clear();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   The spanner is queried for each point, and results stored in the proximity
//   lists.  The query for each point goes out to its cutoff distance, plus the
//   cushion distance.  The isPairDisqualified() function is checked for each
//   pair, to filter out unwanted pairs from the proximity lists.
// Parameters:
//   none
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
template < class _IsPairExcluded, class _PairInteractionCutoff >
void VerletList<_IsPairExcluded, _PairInteractionCutoff>::updateProximityLists() {
  //#pragma omp parallel for
  // For some reason, parallelizing this loop leads to an internal compiler error
  // with no intelligible error message
  for ( size_t id1 = 0; id1 < verletListPoints.size(); id1++ ) {
    double atomCutoffPlusCushion = verletListPoints[id1].getCutoffDist() + cushion;
    vector<int> nearbyPointIDs;
    spanner->getNearbyPoints( id1, atomCutoffPlusCushion, nearbyPointIDs );
    
    vector<VL_NeighborInfo> *neighbors = &proximityLists[id1];
    neighbors->clear();
    
    int externalID1 = verletListPoints[id1].getID();
    int externalID2;
    const Vec3 *vec1 = verletListPoints[id1].getPosition();
    const Vec3 *vec2;
    double cutoff;
    Vec3 delta;
    VL_NeighborInfo neighborInfo;
    int id2;
    double distToCheck;
    for ( size_t j=0; j<nearbyPointIDs.size(); j++ ) {
      id2 = nearbyPointIDs[j];
      externalID2 = verletListPoints[id2].getID();
      if ( isPairExcluded( externalID1, externalID2 ) ) continue;
       
      //further filter results by throwing away pairs that are too far
      cutoff = pairInteractionCutoff( externalID1, externalID2 );
      distToCheck = cutoff + cushion;      
      vec2 = verletListPoints[id2].getPosition();
      delta.x = vec2->x - vec1->x;
      delta.y = vec2->y - vec1->y;
      delta.z = vec2->z - vec1->z;      
      if ( abs(delta.x) > distToCheck || abs(delta.y) > distToCheck || abs(delta.z) > distToCheck ) continue;
      if ( delta.norm2() > distToCheck*distToCheck ) continue;
      
      neighborInfo.id = externalID2;
      neighborInfo.pairCutoff = cutoff;
      neighbors->push_back(neighborInfo);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Update the VerletList.  It will check to see how far points have moved
//   since their last saved positions.  If points have not moved more than
//   cushion/2, then it is done.  The proximity lists are still up-to-date.  But
//   if any point has moved more than cushion/2, then it updates its internal
//   proximity lists, and it stores all the points' current positions as the
//   new saved positions.
// Parameters:
//   none
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
template < class A, class B >
void VerletList<A, B>::update() {
  if (spanner == NULL) {
    cerr << "Error, cannot update an uncommitted VerletList.  Need to call commit() member function." << endl;
    exit(0);
  }

  double movementVerletList = calcMaxMoved();

  // if movement sum exceeds the threshhold,
  // update spanner and then update verlet list
  if ( movementVerletList > movementThreshhold_VerletList ) {
    spanner->update();
    updatePastPositions();
    updateProximityLists();
  }
  // otherwise, if at least the spanner's movement sum
  // exceeds its threshhold, update the spanner
  else {
    double movementSpanner = spanner->calcMaxMovementSinceLastUpdate();
    if ( movementSpanner > movementThreshhold_Spanner ) {
      spanner->update( movementSpanner );
    }
  }
}

template < class A, class B >
void VerletList<A, B>::updatePastPositions() {
  for ( size_t id = 0; id < verletListPoints.size(); id++ ) {
    verletListPoints[id].updatePastPosition();
  }
}

template < class A, class B >
double VerletList<A, B>::calcMaxMoved() const {
  //Compare current point position to the last saved positions.
  //If points have not moved more than cushion/2,
  //then there is no need to update the list
  //of possible contacts
  double maxMoved2 = 0.0;
  for ( size_t id = 0; id < verletListPoints.size(); id++ ) {
    float thisMoved2 = verletListPoints[id].distFromPastPosition2() ;
    if ( thisMoved2 > maxMoved2 ) maxMoved2 = thisMoved2;
  }
  return sqrt( maxMoved2 );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Places the inseted points into a spanner, and establishes needed parameters
// Parameters:
//   beta - spanner parameter beta
//   c - spanner parameter c
//   unitDistance - spanner parameter unitDistance
//   cushion - the VerletList cushion distance
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
template < class A, class B >
void VerletList<A, B>::commit( double beta_, double c_,
                               double interactionCutoff_, double cushion_ ) {
  beta=beta_;
  c=c_;
  cushion = cushion_;
  double spannerBottomEdgeLength = interactionCutoff_ + cushion_; // Angstroms
  double unitDistance = spannerBottomEdgeLength/c; // Angstroms
  delete spanner;
  spanner = new Spanner( beta, c, unitDistance );
  movementThreshhold_VerletList = cushion/2.0;
  movementThreshhold_Spanner = (double)spanner->getMaxAllowedMovement() / 2.0;
  
  for ( size_t id = 0; id < verletListPoints.size(); id++ ) {
    spanner->insert( verletListPoints[id].getPosition(), static_cast<int>(id) );
  }

  //set up proximity lists
  proximityLists.resize( verletListPoints.size() );
  updateProximityLists();
}

template < class A, class B >
void VerletList<A, B>::verify() const {
  spanner->verify();
}

template < class A, class B >
void VerletList<A, B>::printstat() const {
  cout << "Verlet List" << endl;
  cout << "  cushion " << cushion << endl;
  cout << "  Verlet List Movement Threshhold " << movementThreshhold_VerletList << endl;
  cout << "  Spanner Movement Threshhold " << movementThreshhold_Spanner << endl;
  spanner->printStat();
}


#endif 
