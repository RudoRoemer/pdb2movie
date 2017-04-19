#include "VL.h"
#include "RigidUnitSystem.h"

VL::VL( const RigidUnitSystem *rigidUnitSystem,
        const Repulsion &repulsion ) :
  testForExclusion( rigidUnitSystem ),
  verletList( testForExclusion, repulsion )
{
  double maxMaxCutoff = 0.0;
  size_t nPoints = rigidUnitSystem->nPoints();
  for ( size_t p=0; p<nPoints; p++ ) {
    double maxCutoff = repulsion.getMaxInteractionCutoff(p);
    verletList.insert( &rigidUnitSystem->meanPositions(p), maxCutoff, p );
    if ( maxCutoff > maxMaxCutoff ) maxMaxCutoff = maxCutoff;
  }
  double beta = 3.9;
  double c = 4.0;
  double cushion = 2.0;
  verletList.commit( beta, c, maxMaxCutoff, cushion );
  //verletList.printstat();
}

VL::~VL()
{
}

void VL::update() {
  verletList.update();
  notifyObservers();
}

