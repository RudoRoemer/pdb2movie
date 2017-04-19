#ifndef _PROXIMITYMONITOR_H_
#define _PROXIMITYMONITOR_H_

#include "global_defs.h"
#include "MovingPoint.h"
#include "Spanner.h"
#include "SiteID.h"

// class ProximityAtom derives from class MovingPoint.
// It merely adds an interaction cutoff distance parameter,
// used in the ProximityMonitor class.
//////////////////////////////////////////////////////////////////////
class ProximityAtom : public MovingPoint {
public:
  ProximityAtom() : MovingPoint(), cutoffDist(0) {}
  ProximityAtom( const Vector *position, double cutoffDist, SiteID id );
  ProximityAtom( const ProximityAtom& other );
  ~ProximityAtom();
  const ProximityAtom& operator=( const ProximityAtom& other );
  double getCutoffDist() const {return cutoffDist;}
private:
  double cutoffDist;
};

// class ProximityMonitor is used for keeping track of
// pairs of atoms that are within interacting distance.  It uses
// a smart updating feature that only refreshes its internal lists of
// interacting atoms when necessary, described below.  The user interacts
// with the Proximity Monitor as follows:
//   First, insert atoms into the ProximityMonitor, specifying the interaction
//     cutoff distance of each.
//   Commit the ProximityMonitor, which locks the list of atoms and puts the
//     atoms into a spanner.  The user here must supply the desired
//     spanner parameters, and also the desired cushion parameter (see below).
//   Update the proximity monitor as atoms move.
//   Query the ProximityMonitor for the set of atoms that are within the
//     interaction cutoff distance of a particular atom.  You will get some
//     extras here, particularly atoms that were in the "cushion region"
//     (see below).
//
// Here's how the smart updating works.  These details are abstracted away from the user.
// When the user calls the update
// member function, the ProximityMonitor queries the Spanner for a list of
// nearby atoms relative to every atom.  But it does not just query out to the
// interaction cutoff of each atom.  Instead, it queries out to the cutoff plus
// an extra cushion distance.  For each atom, the ProximityMonitor stores a list of nearby atoms,
// which includes any extra atoms that lie in the cushion region.  As long as atoms
// do not move more than cushion/2, this stored list of nearby atoms is guaranteed
// to contain all atoms within the cutoff distance.
//
// When the user subsequently calls update(), the ProximityMonitor will check 
// to see how far atoms have moved.  If atoms have not moved more than cushion/2,
// then it will not even bother querying the spanner.  These update() calls can be
// very fast.  As atoms keep moving, at some point one of the atoms will wander
// more than cutoff/2.  Then, when update() is called, the ProximityMonitor knows that
// its lists of nearby atoms are no longer guaranteed to be accuate, so it re-queries
// the spanner to get up-to-date lists.
//
// This class contains an important virtual function, isPairDisqualified,
// used to derive other classes tailored to specific circumstances.  When the
// ProximityMonitor updates its proximity lists, if you want to exclude particular
// pairs from making it onto the lists, then you write your own isPairDisqualified
// function in a derived class.  For example, if you want to exclude pairs that are the same charge,
// or that are in the same ghost.  This way, you exclude the pair before it makes
// it onto the list.  Otherwise, the pair would be listed in the proximity lists, and you would
// retrieve the pair every time you queried the proximity monitor.  See, for example,
// the derived classes PhobicProximity and PolarProximity.
////////////////////////////////////////////////////////////////////////////////
class ProximityMonitor {

public:
  typedef vector<SiteID> vecID;
  typedef map<SiteID,ProximityAtom> MapIDtoAtom;
  typedef map<SiteID,vecID> ProximityLists;
  
  ProximityMonitor();
  virtual ~ProximityMonitor();

  const vecID& getProximityList( SiteID atom1 ) const;
  void update();
  void free();
  void commit( double beta_, double c_, double unitDistance_, double cushion_ );
  void insert( const Vector *position, double cutoffDist, SiteID id );

protected:
  virtual bool isPairDisqualified( const ProximityAtom &atom1 , const ProximityAtom &atom2 ) const;

private:
  Spanner *spanner;
  ProximityLists proximityLists;
  MapIDtoAtom mapIDtoAtom;
  double cushion;
  float beta;
  float c;
  float unitDistance;
  float triggerUpdateDistance;
  
  void updateProximityLists();
};

////////////////////////////////////////////////////////////////////////////////
// Description:
//   ProximityAtom Constructor.  
// Parameters:
//   position - the address of the position 3-D Vector
//   cutoffDist - the center-to-center distance from the atom beyond which
//     no interactions occur.  
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
inline ProximityAtom::ProximityAtom( const Vector *position_, double cutoffDist_, SiteID id_ ) : 
  MovingPoint( position_, id_ ),
  cutoffDist(cutoffDist_) {}

inline ProximityAtom::ProximityAtom( const ProximityAtom& other ) : MovingPoint(other),
     cutoffDist(other.cutoffDist) {}
////////////////////////////////////////////////////////////////////////////////
// Description:
//   defines the '=' operator for ProximityAtoms.
////////////////////////////////////////////////////////////////////////////////
inline const ProximityAtom& ProximityAtom::operator=(const ProximityAtom& other) {
  position = other.position;
  pastPosition = other.pastPosition;
  id = other.id;
  cutoffDist = other.cutoffDist;
  return other;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   ProximityAtom Destructor.  
////////////////////////////////////////////////////////////////////////////////
inline ProximityAtom::~ProximityAtom() {}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Get the list of atoms that are within the cutoff distance of a particular atom
// Parameters:
//   atom1 - id of particular atom
// Return Value List:
//   returns all atoms that are within the cutoff distance of the input atom.
//   Can also include extra atoms that are in the cushion region, outside the
//   cutoff distance.
////////////////////////////////////////////////////////////////////////////////
inline const ProximityMonitor::vecID& ProximityMonitor::getProximityList( SiteID atom1 ) const {
  ProximityLists::const_iterator iter = proximityLists.find( atom1 );
  if ( iter == proximityLists.end() ) {
    cerr << "Warning, querying contacts of an atom that does not exist " << endl;
    exit(0);
  }
  return iter->second;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   ProximityMonitor Constructor.  
// Parameters:
//   none
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
inline ProximityMonitor::ProximityMonitor() : spanner(NULL),
     cushion(0),
     beta(0),
     c(0),
     unitDistance(0) {
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   ProximityAtom Destructor.  
////////////////////////////////////////////////////////////////////////////////
inline ProximityMonitor::~ProximityMonitor() {
  delete spanner;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Insert an atom into the ProximityMonitor.  Atoms can be inserted only when
//   setting up the ProximityMonitor.  After you call commit(), you cannot
//   insert atoms anymore.
// Parameters:
//   position - the address of the point's 3-D Vector.  As the point moves,
//     the proximity monitor will always know what its new position is.
//     IMPORTANT - this address must always be valid!  If the 3-D Vector
//     ends up moving to a different location in memory, this will be bad.
//     For example, if the 3-D Vector is stored in a std::vector<>, and then
//     you extend the size of the std::vector, the std::vector may move its
//     contents somewhere else, invalidating all pointers to its contents.
//   cutoffDist - defines the center-to-center distance beyond which atoms are
//     not neighbors.  When you query the ProximityMonitor for this atom, it will return atoms
//     that are within cutoffDist of this atom.
//   id - the ID of the point
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
inline void ProximityMonitor::insert( const Vector *position, double cutoffDist, SiteID id ) {

  if ( spanner != NULL ) {
    cerr << "Error, cannot insert atom into committed proximity monitor" << endl;
    exit(0);
  }
  
  MapIDtoAtom::const_iterator iter = mapIDtoAtom.find( id );
  if ( iter != mapIDtoAtom.end() ) {
    cerr << "Error, attempting to insert duplicate atom into proximity monitor " << endl;
    exit(0);
  }
  
  mapIDtoAtom[id] = ProximityAtom( position, cutoffDist, id );

  //set up a blank entry in the proximity lists
  proximityLists[id];
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Clears the internal spanner, the stored atoms, and the stored
//   proximity lists
// Parameters:
//   none
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
inline void ProximityMonitor::free() {
  delete spanner;
  spanner = NULL;
  mapIDtoAtom.clear();
  proximityLists.clear();
}

#endif 
