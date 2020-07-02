#include <cmath>
#include <iostream>

#include "Util.h"

namespace Util {

double alphas(double t, bool dla, double q0) {
  if (dla) {
    return 0.3;
  }
  else {
    double beta0 = ( 11.0 * Nc - 2.0 * Nf ) / 3.0;
    if ( t > q0 ) return std::min( 2.0 * M_PI / beta0 / log( t / q0 ), 1.0 );
    else return 1.0;
  }
}

double Pgg(double z, bool dla) {
  if ( dla ) return CA / z;
  return CA * pow( 1.0 - z * (1.0 - z), 2.0 ) / z / (1.0 - z);
}

double Pgq(double z, bool dla) {
  if ( dla ) return CF / z;
  else return CF * ( 1.0 + (1.0 - z) * (1.0 - z) ) / z; 
}

double Pqq(double z, bool dla) {
  return CF * ( 1.0 + z * z ) / ( 1.0 - z ); 
}

double Pqg(double z, bool dla) {
  return Nf / 2.0 * ( z * z + ( 1.0 - z ) * ( 1.0 - z ) ); 
}

} // end namespace Util
