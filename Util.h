#ifndef Util_H
#define Util_H

#ifndef hbarc
#define hbarc (0.197327)
#endif

#ifndef egamma
#define egamma (0.5772156649)
#endif

#ifndef CA
#define CA (3.0)
#endif

#ifndef CF
#define CF (4.0/3.0)
#endif

#ifndef Nc
#define Nc (3.0)
#endif

#ifndef Nf
#define Nf (5.0)
#endif

#ifndef eps1
#define eps1 (0.001)
#endif

#include <string>

namespace Util {

  std::string StringFind4(std::string file_name, std::string str_in);

  // Splitting functions
  double Pgg(double z, bool dla);
  double Pgq(double z, bool dla);
  double Pqq(double z, bool dla);
  double Pqg(double z, bool dla);
    
  // Running alpha_s
  double alphas(double t, bool dla, double q0);

}

#endif // Util_H
