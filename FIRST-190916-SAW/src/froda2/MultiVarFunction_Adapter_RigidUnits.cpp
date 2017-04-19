#include "MultiVarFunction_Adapter_RigidUnits.h"
#include "RigidUnitSystem.h"
#include "GeneralizedCoords.h"
#include "ConstraintEnforcingPotential.h"
#include <cmath>

MultiVarFunction_Adapter_RigidUnits::
  MultiVarFunction_Adapter_RigidUnits(
    RigidUnitSystem *rigidUnitSystem_,
    ConstraintEnforcingPotential *cep_ ) :
  rigidUnitSystem( rigidUnitSystem_ ),
  cep( cep_ )
{
}

MultiVarFunction_Adapter_RigidUnits::
  MultiVarFunction_Adapter_RigidUnits() :
  rigidUnitSystem( NULL ),
  cep( NULL )
{
}
  
MultiVarFunction_Adapter_RigidUnits::~MultiVarFunction_Adapter_RigidUnits()
{
}

void MultiVarFunction_Adapter_RigidUnits::
  setRigidUnits(
    RigidUnitSystem *rigidUnitSystem_,
    ConstraintEnforcingPotential *cep_ )
{
  rigidUnitSystem = rigidUnitSystem_;
  cep = cep_;
}

void MultiVarFunction_Adapter_RigidUnits::getQ( std::vector<double> &q )
{
  int nRU = rigidUnitSystem->nRigidUnits();
  q.resize( nRU*6 );
  int i = 0;
  for ( int ru=0;  ru < nRU; i += 3, ru++ ) {
    q[i] = rigidUnitSystem->centers(ru).x;
    q[i+1] = rigidUnitSystem->centers(ru).y;
    q[i+2] = rigidUnitSystem->centers(ru).z;
  }
  for ( int ru=0;  ru < nRU; i += 3, ru++ ) {
    //here we scale the rotors
    q[i] = rigidUnitSystem->rotors(ru).x*rigidUnitSystem->radius(ru);
    q[i+1] = rigidUnitSystem->rotors(ru).y*rigidUnitSystem->radius(ru);
    q[i+2] = rigidUnitSystem->rotors(ru).z*rigidUnitSystem->radius(ru);
  }
}

void MultiVarFunction_Adapter_RigidUnits::setQ( const std::vector<double> &q )
{
  size_t nRU = rigidUnitSystem->nRigidUnits();
  if ( q.size() != nRU*6 ) {
    cout << "Error: Q must have size nRigidUnits*6 " << endl;
    exit(0);
  }
  size_t i = 0; //counts over the centers 
  size_t j = 3*nRU;
  for ( size_t ru=0;  ru < nRU; i += 3, j += 3, ru++ ) {
    if ( rigidUnitSystem->hasZeroRadius(ru) ) {
      rigidUnitSystem->setCenter( ru, 
                         Vec3( q[i], q[i+1], q[i+2] ) );
    }
    else {
      rigidUnitSystem->setCenterAndRotor( ru, 
                                 Vec3( q[i], q[i+1], q[i+2] ),
                                 Vec3( q[j]/rigidUnitSystem->radius(ru), 
                                       q[j+1]/rigidUnitSystem->radius(ru),
                                       q[j+2]/rigidUnitSystem->radius(ru) ) );
    }
  }
  rigidUnitSystem->update();
}

double MultiVarFunction_Adapter_RigidUnits::eval()
{
  return cep->energy();
}

void MultiVarFunction_Adapter_RigidUnits::getNegGrad( std::vector<double> &r )
{
  const GeneralizedCoords *grad = &cep->gradient();
  size_t nRU = rigidUnitSystem->nRigidUnits();
  r.resize( nRU*6 );
  size_t i = 0;
  for ( size_t ru=0;  ru < nRU; i += 3, ru++ ) {
    r[i] = -grad->centersRU[ru].x;
    r[i+1] = -grad->centersRU[ru].y;
    r[i+2] = -grad->centersRU[ru].z;
  }
  for ( size_t ru=0;  ru < nRU; i += 3, ru++ ) {
    if ( rigidUnitSystem->hasZeroRadius(ru) ) {
      r[i] = r[i+1] = r[i+2] = 0.0;
    }
    else {
      r[i] = -grad->rotorsRU[ru].x/rigidUnitSystem->radius(ru);
      r[i+1] = -grad->rotorsRU[ru].y/rigidUnitSystem->radius(ru);
      r[i+2] = -grad->rotorsRU[ru].z/rigidUnitSystem->radius(ru);
    }
  }

}

void MultiVarFunction_Adapter_RigidUnits::getPreconditionedNegGrad( std::vector<double> &s, bool &isPreconditioned )
{
  //In this function, the preconditioner (d2V_dQ2) is recalculated.
  //So, every time getPreconditionedNegGrad is called during
  //the minimization, the preconditioner itself changes
  
  getNegGrad( s );

  const GeneralizedCoords *d2V_dQ2_diagonal = 
    &cep->d2V_dQ2_diagonal();
  
  size_t nRU = rigidUnitSystem->nRigidUnits();
  for ( size_t ru=0; ru<nRU; ru++ ) {
    if ( d2V_dQ2_diagonal->centersRU[ru].x < numeric_limits<double>::epsilon() ||
         d2V_dQ2_diagonal->centersRU[ru].y < numeric_limits<double>::epsilon() ||
         d2V_dQ2_diagonal->centersRU[ru].z < numeric_limits<double>::epsilon() ||
         d2V_dQ2_diagonal->rotorsRU[ru].x < -numeric_limits<double>::epsilon() ||
         d2V_dQ2_diagonal->rotorsRU[ru].y < -numeric_limits<double>::epsilon() ||
         d2V_dQ2_diagonal->rotorsRU[ru].z < -numeric_limits<double>::epsilon() )
    {
      isPreconditioned = false;
      return;
    }
  }
  
  double rotorScaleFactor;
  size_t i = 0;
  for ( size_t ru=0;  ru < nRU; i += 3, ru++ ) {
    s[i] /= d2V_dQ2_diagonal->centersRU[ru].x;
    s[i+1] /= d2V_dQ2_diagonal->centersRU[ru].y;
    s[i+2] /= d2V_dQ2_diagonal->centersRU[ru].z;
  }
  for ( size_t ru=0;  ru < nRU; i += 3, ru++ ) {
    rotorScaleFactor = rigidUnitSystem->radius(ru)*rigidUnitSystem->radius(ru);
    s[i] = s[i]/(d2V_dQ2_diagonal->rotorsRU[ru].x/rotorScaleFactor + 1.0);
    s[i+1] = s[i+1]/(d2V_dQ2_diagonal->rotorsRU[ru].y/rotorScaleFactor + 1.0);
    s[i+2] = s[i+2]/(d2V_dQ2_diagonal->rotorsRU[ru].z/rotorScaleFactor + 1.0);
  }  
  isPreconditioned = true;
}

