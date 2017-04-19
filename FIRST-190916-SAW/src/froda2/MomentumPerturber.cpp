#include "MomentumPerturber.h"
#include "RigidUnitSystem.h"

MomentumPerturber::MomentumPerturber( RigidUnitSystem *rigidUnitSystem_ ) :
  rigidUnitSystem( rigidUnitSystem_ ),
  rigidUnitFitter( rigidUnitSystem ),
  isQ1set( false ),
  readyToPerturb( false )
{  
}

MomentumPerturber::~MomentumPerturber()
{
}

void MomentumPerturber::setQ1() { 
  Q1_rigidUnitPoints = rigidUnitSystem->absolutePositions();
  isQ1set = true;  
}

void MomentumPerturber::determineDeltaQ() {
  if ( !isQ1set ) return;
  
  rigidUnitSystem->collapseRotors();
  rigidUnitFitter.calcFitToRigidUnitPoints( Q1_rigidUnitPoints );
  
  isQ1set = false;
  readyToPerturb = true;
}

void MomentumPerturber::perturb() {
  if ( !readyToPerturb ) return;
  
  size_t nRU = rigidUnitSystem->nRigidUnits();
  const double maglim = 0.5;
  const double maglim2 = maglim*maglim;
  double mag2;
  
  rigidUnitSystem->collapseRotors();

  for ( size_t ru = 0; ru < nRU; ru++ ) {
    //get the net translation of this rigid unit from state Q1 to state Q2
    Vec3 deltaCenter = rigidUnitFitter.getFitTranslation(ru);
    //we need to set it equal to its negative, because the fit
    //was from the state Q2 onto Q1, instead of Q1 onto Q2
    deltaCenter.x = -deltaCenter.x; 
    deltaCenter.y = -deltaCenter.y; 
    deltaCenter.z = -deltaCenter.z; 

    //if translation was too large, truncate
    mag2 = deltaCenter.norm2();
    if ( mag2 > maglim2 ) {
      deltaCenter *= maglim/sqrt(mag2);
    }

    //apply perturbation
    rigidUnitSystem->addToCenter( ru, deltaCenter );
    
    if ( rigidUnitSystem->hasZeroRadius(ru) ) continue;
    
    //get the net rotation of this rigid unit from state Q1 to Q2
    Vec3 rotor = rigidUnitFitter.getFitRotation(ru);
    rotor.x = -rotor.x;
    rotor.y = -rotor.y;
    rotor.z = -rotor.z;
    
    //if rotation was too large, truncate.
    //use the rotor magnitude (which is the rotation angle for small angles),
    //scaled by the radius
    mag2 = rotor.norm2()*rigidUnitSystem->radius(ru)*rigidUnitSystem->radius(ru);
    if ( mag2 > maglim2 ) {
      rotor *= maglim/sqrt(mag2);
    }
    
    // apply perturbation
    rigidUnitSystem->setRotor( ru, rotor );
  }
  rigidUnitSystem->collapseRotors();
  rigidUnitSystem->update();
  
  readyToPerturb = false;
}
