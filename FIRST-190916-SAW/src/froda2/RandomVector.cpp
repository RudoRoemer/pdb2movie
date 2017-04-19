#include "RandomVector.h"
#include "mt19937ar.h"
#include "Vect.h"
typedef Vect<double> Vec3;
#include <cmath>

const double PI = 3.141592653589793238462643383279;
const double TWOPI(2.0 * PI);

void generateRandomUnitVector( Vec3& vec ) {
  double phi = genrand_real2()*TWOPI;
  //double theta = genrand_real2()*PI;
  //double costheta = cos(theta);
  //double sintheta = sin(theta);
  double costheta = genrand_real1()*2.0 - 1.0;
  double sintheta = sqrt(1.0 - costheta*costheta);
  vec.x = sintheta*cos(phi);
  vec.y = sintheta*sin(phi);
  vec.z = costheta;
}
