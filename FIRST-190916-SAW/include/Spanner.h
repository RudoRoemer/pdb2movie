// Created by Dan Farrell
//   June 16, 2006

#ifndef SPANNER_H_
#define SPANNER_H_

#include "DiscreteHierarchy.h"
#include "MovingPoint.h"
#include "Vect.h"
typedef Vect<double> Vec3;
#include <limits>
//#include "Timer.h"

typedef DiscreteHierarchy<MovingPoint> DCH;

// The Spanner class provides an interface for using a Discrete Centers Heirarchy (dch).
// The user of the Spanner class can do the following:
//   insert point positions along with corresponding point ids,
//   update the Spanner periodically as the points move,
//   ask the Spanner for a list of near neighbors of any of the points.
// Behind the scenes, the Spanner watches how far points move and translates that movement into 
// a "time" that is used by the dch.  The Spanner also checks to see if points move
// more than the max allowed distance.  If this happens, the Spanner resets the dch.
// For more information, see the individual member function definitions.
class Spanner {
public:
  typedef map<int, MovingPoint> MapIDtoPoint;
  
  Spanner( float beta, float c, float unitDistance );
  ~Spanner();
  void insert( const Vec3 *position, int id );
  void remove( int id );
  void rebuild();
  void update();
  void update( float maxMoved );
  float calcMaxMovementSinceLastUpdate() const;
  void printStat() const;
  void verify() const;
  void getNearbyPoints( int id, float targetDist, vector<int>& nearbyPointIDs ) const;
  float getBeta() const {return beta;}
  float getC() const {return c;}
  float getUnitDistance() const {return unitDistance;}
  float getMaxAllowedMovement() const {return maxAllowedMovement;}
  static float getMaxAllowedMovement( float beta, float c, float unitDistance );
private:
  MapIDtoPoint mapIDtoPoint;
  float time;
  float beta;
  float c;
  float unitDistance;
  float timeResolution;
  float maxTime;
  DCH *dch;
  float maxAllowedMovement;
  //Timer timer;
  vector<const MovingPoint*> nearbyPoints;
  vector<int> nearbyPointIDs;
};

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Finds points that are closer than a given distance to the input point.
//   It does this by querying the dch bottom-level neighbors.
// Parameters:
//   id - the id of the point whose neighbors you are requesting
//   targetDist - points that are closer than this distance will be output,
//      with some extra points that are actually farther than this distance
// Return Value List:
//   a vector containing the neighbors of p
////////////////////////////////////////////////////////////////////////////////
inline void Spanner::getNearbyPoints( 
    int id, 
    float targetDist, 
    vector<int>& nearbyPointIDs ) const {
  
  MapIDtoPoint::const_iterator iter = mapIDtoPoint.find( id );
  if ( iter == mapIDtoPoint.end() ) {
    cerr << "Warning, querying a point that does not exist in spanner " << endl;
    exit(0);
  }
    
  const MovingPoint *point = &iter->second;
  dch->getBottomLevelNeighbors( point, targetDist, nearbyPointIDs );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Output various statistics of the dch to cout
// Parameters:
//   none
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
inline void Spanner::printStat() const {
  dch->printstat();
  cout << "Spanner: " << endl;
  cout << "  Max allowed movement is " << getMaxAllowedMovement() << endl;
  cout << "  Certificate time resolution is " << timeResolution << endl;
  cout << "  Max Time is  " << maxTime << endl;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Constructor.  Takes the input dch parameters and instantiates a
//   Discrete Centers Heirarchy (dch).
// Parameters:
//   (See dch documentation for info on these dch parameters)
//   beta_ - dch parameter beta
//   c_ - dch parameter c
//   unitDistance_ - dch parameter minRadius
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
inline Spanner::Spanner( float beta_, float c_, float unitDistance_) : time(0),
                                                                   beta(beta_),
                                                                   c(c_),
                                                                   unitDistance(unitDistance_),
                                                                   dch(NULL) {
  maxAllowedMovement = getMaxAllowedMovement( beta, c, unitDistance );
  //timeResolution = 0x1p-8; //hex for 1 * 2^(-8), which is ~ 0.003
  //maxTime = (0x1p24 - 0x1p0) * timeResolution; // to keep the desired resolution, we must not let time exceed this value
  timeResolution = 0.0001;
  maxTime = timeResolution/numeric_limits<float>::epsilon();
  dch = new DCH( beta, c, unitDistance );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Destructor.
// Parameters:
//   none
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
inline Spanner::~Spanner() {
  delete dch;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Insert a point into the spanner.  Points can be inserted or removed at
//   any time, but the spanner must first be up-to-date.  Behavior is not
//   defined if you insert a point into a non-updated spanner.
// Parameters:
//   position - the address of the point's 3-D Vector.  As the point moves,
//     the spanner will always know what its new position is.
//     IMPORTANT - this address must always be valid!  If the 3-D Vector
//     ends up moving to a different location in memory, this will be bad.
//     For example, if the 3-D Vector is stored in a std::vector<>, and then
//     you extend the size of the std::vector, the std::vector may move its
//     contents somewhere else, invalidating all pointers to its contents.
//   id - the ID of the point
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
inline void Spanner::insert( const Vec3 *position, int id ) {
  MapIDtoPoint::const_iterator iter = mapIDtoPoint.find( id );
  if ( iter != mapIDtoPoint.end() ) {
    cerr << "Error, attempting to insert duplicate point into spanner " << endl;
    exit(0);
  }
  mapIDtoPoint[id] = MovingPoint( position, id );
  dch->insert( &mapIDtoPoint[id] );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Remove a point from the spanner.  Points can be inserted or removed at
//   any time, but the spanner must first be up-to-date.  Behavior is not
//   defined if you remove a point from a non-updated spanner.
// Parameters:
//   id - the id of the point you wish to remove
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
inline void Spanner::remove( int id ) {
  MapIDtoPoint::iterator iter = mapIDtoPoint.find( id );
  if ( iter == mapIDtoPoint.end() ) {
    cerr << "Error, attempting to remove point from spanner that does not exist" << endl;
    exit(0);
  }
  dch->remove( &iter->second );
  mapIDtoPoint.erase( iter );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Compute the maximum allowed movement between updates.  This function is
//   general, meaning that you can compute the max allowed movement
//   for any spanner parameters you want.  You don't have to use the
//   parameters of this particular spanner object.
// Parameters:
//   beta - dch parameter beta
//   c - dch parameter c
//   unitDistance - dch parameter minRadius
// Return Value List:
//   the maximum allowed movement between updates
////////////////////////////////////////////////////////////////////////////////
inline float Spanner::getMaxAllowedMovement( float beta, float c, float unitDistance ) {
  return (c*(beta - 1.0f) - 2.0f*beta)/4.0f*unitDistance;
}
#endif /*SPANNER_H_*/
