// Created by Dan Farrell
//   Jun 16, 2006

#include "Spanner.h"
#include "MovingPoint.h"
#include <cmath>
#include <iomanip>
#include <iostream>

// initialize the static member variable idcount to 0.
template <class MP>
int Node<MP>::idcount = 0;

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Delete the current Discrete Centers Heirarchy and rebuild it.  You do not
//   need to reinsert the points.  The Spanner object keeps the points you
//   have inserted, and re-inserts them for you into a fresh, clean dch.
// Parameters:
//   id - the id of the point you wish to remove
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
void Spanner::rebuild() {
  //reset the time to 0.
  time = 0;
  
  //delete and rebuild the dch
  delete dch;
  Node<MovingPoint>::resetIDcount();
  dch = new DCH( beta, c, unitDistance );
  for ( MapIDtoPoint::iterator iter = mapIDtoPoint.begin(); iter != mapIDtoPoint.end(); iter++ ) {
    MovingPoint *point = &iter->second;
    point->updatePastPosition();
    dch->insert( point );
  }
}

void Spanner::verify() const {
  for ( MapIDtoPoint::const_iterator iter = mapIDtoPoint.begin(); iter != mapIDtoPoint.end(); iter++ ) {
    const MovingPoint *point = &iter->second;
    if ( point->getPosition() == NULL) {
      cout << "Error: spanner point's position is null" << endl;
      exit(0);
    }
    if ( isnan(point->getPosition()->x) || 
         isnan(point->getPosition()->y) ||
         isnan(point->getPosition()->z) ) {
      cout << "Found Nan: Spanner point " << point->getID() 
           << " " << *point->getPosition() << endl;
      exit(0);
    }
  }
  dch->verify();
}
  

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Update the spanner and dch.  Note that there are multiple overloaded versions of
//   this function.  This one takes no inputs.  It first looks at how far
//   the points have moved, and then translates that into a "time"
//   required by the dch.
// Parameters:
//   none
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
void Spanner::update() {
  //cout << "Spanner: maxMoved = " << maxMoved << endl;
  update( calcMaxMovementSinceLastUpdate() );
}

float Spanner::calcMaxMovementSinceLastUpdate() const {
  //Compare current point position to the last saved positions.
  float maxMoved2=0;
  for ( MapIDtoPoint::const_iterator iter = mapIDtoPoint.begin(); iter != mapIDtoPoint.end(); iter++ ) {
    const MovingPoint *point = &iter->second;
    float thisMoved2 =  point->distFromPastPosition2() ;
    if ( thisMoved2 > maxMoved2 ) maxMoved2 = thisMoved2;
  }
  return sqrt( maxMoved2 );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Update the spanner and dch.  Note that there are multiple overloaded versions of
//   this function.  In this one, you pass in the maximum distance that any point
//   has moved since the last update.  This is helpful in cases where you
//   already know the max distance points have moved, and you have no need to 
//   have the spanner recompute it.
// Parameters:
//   maxMoved - the maximum distance that any point in the spanner has moved
//   since the last spanner update.
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
void Spanner::update( float maxMoved ) {
  if ( maxMoved > maxAllowedMovement ) {
    //cout << "  Spanner max allowed movement was exceeded.  Regenerating spanner. " << endl;
    rebuild();
    //cout << "Done regenerating spanner." << endl;
    return;
  }

  time += maxMoved;
  if (time > maxTime ){
    //std::cout << "Spanner internal ticker has exceed max value, resetting spanner" << std::endl;
    rebuild();
    //cout << "Done regenerating spanner." << endl;
  }
  for ( MapIDtoPoint::iterator iter = mapIDtoPoint.begin(); iter != mapIDtoPoint.end(); iter++ ) {
    MovingPoint *point = &iter->second;
    point->updatePastPosition();
  }
  dch->update( time );
  //cout << setprecision(30) << "Verifying dch after update of " << maxMoved << endl;
  //cout << setprecision(30) << "Spanner Time is " << time << endl;
  //cout << setprecision(30) << "dch Time is " << dch->getTime() << endl;
  //dch->verify();
}


