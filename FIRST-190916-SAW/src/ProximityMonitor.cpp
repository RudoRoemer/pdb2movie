#include "global_defs.h"
#include "ProximityMonitor.h"
#include "MolFramework.h"
#include "Spanner.h"

extern Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
// Description:
//   When the ProximityMonitor is updating its proximity lists, this function
//   is checked to see if certain atom-atom pairs should be disqualified
//   from being added to the list.  In this base class, it returns false,
//   meaning that no pair will be disqualified.  But in derived classes,
//   you can make it check any criteria you want, such as whether the
//   atoms are 3rd neighbors, or in the same ghost, or have like charges.
//   By the way, when I say derived class, I don't mean make a copy of this class
//   modify it!  I mean make a new class that inherits from this class.
// Parameters:
//   atom1, atom2 - the pair of atoms under consideration
// Return Value List:
//   true - the pair is to be disqualified
//   false - the pair is not to be disqualified
////////////////////////////////////////////////////////////////////////////////
bool ProximityMonitor::isPairDisqualified( const ProximityAtom &atom1 , const ProximityAtom &atom2 ) const {
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   The spanner is queried for each atom, and results stored in the proximity
//   lists.  The query for each atom goes out to its cutoff distance, plus the
//   cushion distance.  The isPairDisqualified() function is checked for each
//   pair, to filter out unwanted pairs from the proximity lists.
// Parameters:
//   none
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
void ProximityMonitor::updateProximityLists() {
  
  vector<int> nearbyAtomIDs;
  Vector delta;
  
  for ( MapIDtoAtom::const_iterator iter = mapIDtoAtom.begin(); iter != mapIDtoAtom.end(); iter++ ) {
    const ProximityAtom *atom1 = &iter->second;
    const ProximityAtom *atom2;
    SiteID id1 = atom1->getID();
    vecID *neighbors = &proximityLists[id1];
    neighbors->clear();
    double distToCheck = atom1->getCutoffDist() + cushion;
    spanner->getNearbyPoints( id1, distToCheck, nearbyAtomIDs );     
    
    for ( size_t j=0; j<nearbyAtomIDs.size(); j++ ) {
      SiteID id2 = nearbyAtomIDs[j];
      atom2 = &mapIDtoAtom[ id2 ];
       
      //further filter results by throwing away pairs that are too far
      delta = *atom1->getPosition() - *atom2->getPosition();
      if ( abs(delta.x) > distToCheck || abs(delta.y) > distToCheck || abs(delta.z) > distToCheck ) continue;
      if ( isPairDisqualified( *atom1, *atom2 ) ) continue;

      neighbors->push_back(id2);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Update the ProximityMonitor.  It will check to see how far atoms have moved
//   since their last saved positions.  If atoms have not moved more than
//   cushion/2, then it is done.  The proximity lists are still up-to-date.  But
//   if any atom has moved more than cushion/2, then it updates its internal
//   proximity lists, and it stores all the atoms' current positions as the
//   new saved positions.
// Parameters:
//   none
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
void ProximityMonitor::update() {
  if (spanner == NULL) {
    cerr << "Error, cannot update an uncommitted ProximityMonitor.  Need to call commit() member function." << endl;
    exit(0);
  }
  //Compare current atom position to the last saved positions.
  //If atoms have not moved more than cushion/2,
  //then there is no need to update the list
  //of possible contacts
  double maxMoved2=0;
  for ( MapIDtoAtom::const_iterator iter = mapIDtoAtom.begin(); iter != mapIDtoAtom.end(); iter++ ) {
    const ProximityAtom *atom = &iter->second;
    float thisMoved2 =  atom->distFromPastPosition2() ;
    if ( thisMoved2 > maxMoved2 ) maxMoved2 = thisMoved2;
  }
  double maxMoved = sqrt( maxMoved2 );
  if ( maxMoved > triggerUpdateDistance ) {
    spanner->update( maxMoved );
    for ( MapIDtoAtom::iterator iter = mapIDtoAtom.begin(); iter != mapIDtoAtom.end(); iter++ ) {
      ProximityAtom *atom = &iter->second;
      atom->updatePastPosition();
    }
    updateProximityLists();
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Places the inseted atoms into a spanner, and establishes needed parameters
// Parameters:
//   beta - spanner parameter beta
//   c - spanner parameter c
//   unitDistance - spanner parameter unitDistance
//   cushion - the ProximityMonitor cushion distance
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
void ProximityMonitor::commit( double beta_, double c_, double unitDistance_, double cushion_ ) {
  beta=beta_;
  c=c_;
  unitDistance=unitDistance_;
  cushion = cushion_;
  delete spanner;
  spanner = new Spanner( beta, c, unitDistance );
  triggerUpdateDistance = cushion/2.0;
  if ( triggerUpdateDistance > (double)spanner->getMaxAllowedMovement() ) {
    cerr << "Warning, the amount of movement that triggers a spanner update " << endl;
    cerr << "exceeds the maximum movement allowed by the given spanner parameters." << endl;
  }
  
  for ( MapIDtoAtom::const_iterator iter = mapIDtoAtom.begin(); iter != mapIDtoAtom.end(); iter++ ) {
    const ProximityAtom *atom = &iter->second;
    spanner->insert( atom->getPosition(), atom->getID() );
  }
  if ( parameters.verbose >= 2 ) {
    spanner->printStat( );
  }
  updateProximityLists();
}
