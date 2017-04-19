#include "RigidUnitFitter.h"
#include <iostream>
#include "RigidUnitSystem.h"

using namespace std;


RigidUnitFitter::RigidUnitFitter( const RigidUnitSystem *rigidUnitSystem_ ) :
  rigidUnitSystem( rigidUnitSystem_ )
{
}

RigidUnitFitter::~RigidUnitFitter()
{
}

void RigidUnitFitter::setRigidUnits( const RigidUnitSystem *rigidUnitSystem_ ) {
  rigidUnitSystem = rigidUnitSystem_;
}

void RigidUnitFitter::calcFitToMeanPoints( const vector<Vec3> &target ) {
  if ( target.size() != rigidUnitSystem->nPoints() ) {
    cout << "Error: target/source sizes do not match" << endl;
    exit(0);
  }
  
  size_t nRU = rigidUnitSystem->nRigidUnits();
  vector <Vec3> ruBase;
  vector <Vec3> ruTarget;
  Fit fit;
  fitRotation.resize( nRU );
  fitTranslation.resize( nRU );
  const vector<int> *ruplist;
  size_t nRUP;
  int rup;
  int p;
  for ( size_t ru = 0; ru < nRU; ru++ ) {
    ruplist = &rigidUnitSystem->getRUPlistFromRU( ru );
    nRUP = ruplist->size();
    ruBase.resize( nRUP );
    ruTarget.resize( nRUP );
    for ( size_t i = 0; i < nRUP; i++ ) {
      rup = (*ruplist)[i];
      ruBase[i] = rigidUnitSystem->basePositions(rup);
      p = rigidUnitSystem->getPfromRUP( rup );
      ruTarget[i] = target[p];
    }
    
    fit.setSourceBasePoints( ruBase, rigidUnitSystem->radius(ru) );
    fit.setTargetAbsolutePoints( ruTarget );
    fit.simpleFit();
    fitRotation[ru] = fit.getFitRotor();
    fitTranslation[ru] = fit.getFitCenter() - rigidUnitSystem->centers(ru);
  }
}

void RigidUnitFitter::calcFitToRigidUnitPoints( const vector<Vec3> &absolutePositions ) {
  if ( absolutePositions.size() != rigidUnitSystem->nRigidUnitPoints() ) {
    cout << "Error: target/source sizes do not match" << endl;
    exit(0);
  }

  size_t nRU = rigidUnitSystem->nRigidUnits();
  fitRotation.resize( nRU );
  fitTranslation.resize( nRU );
  #pragma omp parallel
  {
    vector <Vec3> ruBase;
    vector <Vec3> ruTarget;
    Fit fit;;
    const vector< int > *ruplist;
    size_t nRUP;
    int rup;
    #pragma omp for
    for ( size_t ru = 0; ru < nRU; ru++ ) {
      ruplist = &rigidUnitSystem->getRUPlistFromRU( ru );
      nRUP = ruplist->size();
      ruBase.resize( nRUP );
      ruTarget.resize( nRUP );
      for ( size_t i = 0; i < nRUP; i++ ) {
        rup = (*ruplist)[i];
        ruBase[i] = rigidUnitSystem->basePositions(rup);
        ruTarget[i] = absolutePositions[rup];
      }
      
      fit.setSourceBasePoints( ruBase, rigidUnitSystem->radius(ru) );
      fit.setTargetAbsolutePoints( ruTarget );
      fit.simpleFit();
      fitRotation[ru] = fit.getFitRotor();
      fitTranslation[ru] = fit.getFitCenter() - rigidUnitSystem->centers(ru);
    }
  } // end parallel region
}
