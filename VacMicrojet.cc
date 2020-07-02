#include <vector>

#include "VacMicrojet.h"

VacMicrojet::VacMicrojet(double Rjet)
{

  _moment = 1.0;
  _spectrum.clear();
  _Rjet = Rjet;

}

void VacMicrojet::add_spectrum_bin(double Ejet, double moment)
{

  _spectrum.push_back( std::make_pair( Ejet, moment ) );

}
