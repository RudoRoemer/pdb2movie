#include "MultiVarFunction_Adapter_Fit.h"
#include "Fit.h"

MultiVarFunction_Adapter_Fit::MultiVarFunction_Adapter_Fit( Fit *fit_ ) :
  fit( fit_ )
{
}

MultiVarFunction_Adapter_Fit::~MultiVarFunction_Adapter_Fit()
{
}

void MultiVarFunction_Adapter_Fit::getQ( std::vector<double> &q ) {
  Vec3 rotor;
  rotor = fit->getRotor();
  q.resize(3);
  q[0] = rotor.x * fit->getRadiusGyration();
  q[1] = rotor.y * fit->getRadiusGyration();
  q[2] = rotor.z * fit->getRadiusGyration(); 
}

void MultiVarFunction_Adapter_Fit::setQ( const std::vector<double> &q ) {
  Vec3 rotor( q[0]/fit->getRadiusGyration(),
                       q[1]/fit->getRadiusGyration(),
                       q[2]/fit->getRadiusGyration() );
  fit->setRotor( rotor );
}

double MultiVarFunction_Adapter_Fit::eval() {
  return fit->energy();
}

void MultiVarFunction_Adapter_Fit::getNegGrad( std::vector<double> &r ) {
  const Vec3 *grad = &fit->gradient();
  r.resize(3);
  r[0] = -grad->x;
  r[1] = -grad->y;
  r[2] = -grad->z; 
  //cout << grad << endl;
}


void MultiVarFunction_Adapter_Fit::getPreconditionedNegGrad( std::vector<double> &s, bool &isPreconditioned ) {
  const Vec3 *grad = &fit->gradient();
  s.resize(3);
  s[0] = -grad->x;
  s[1] = -grad->y;
  s[2] = -grad->z; 
  
  const Vec3 *d2V_dB2 = &fit->secondDerivativeDiagonal();

  bool isPositiveDefinite = !(
         d2V_dB2->x < numeric_limits<double>::epsilon() ||
         d2V_dB2->y < numeric_limits<double>::epsilon() ||
         d2V_dB2->z < numeric_limits<double>::epsilon() );

  if ( isPositiveDefinite ) {
    s[0] /= d2V_dB2->x;
    s[1] /= d2V_dB2->y;
    s[2] /= d2V_dB2->z;
  }
  isPreconditioned = isPositiveDefinite;
  
}

