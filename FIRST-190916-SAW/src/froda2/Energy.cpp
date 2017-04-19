#include "Energy.h"
#include <list>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "RigidUnitSystem.h"

using namespace std;

Energy::Energy() : 
  isUpdated( false ),
  Etot(0.0)
{
}

Energy::~Energy() {
}

void Energy::receiveNotification( Observable *observable )
{
  const RigidUnitSystem *rigidUnitSystem = 
    dynamic_cast<const RigidUnitSystem *>( observable );
  if ( rigidUnitSystem == NULL ) {
    cout << "Error: dynamic cast failed" << endl;
    exit(0);
  }
  if ( rigidUnitSystem->AbsolutePositionsChanged() ) isUpdated = false;
}

void Energy::update() {
  Etot = 0.0;
  for ( std::list<EnergyTerm*>::iterator it = energyTerms.begin();
        it != energyTerms.end();
        it++ ) {
    // 'it' is an iterator to a pointer
    // *it is a pointer
    Etot += (*it)->energy();  
  }
  /*for (size_t i = 0; i < energyTerms.size(); i++) {
    Etot += (*energyTerms[i])->energy();
    }*/
  isUpdated = true;
  //cout << Etot << " " << syncMovementCounter << endl;
}

