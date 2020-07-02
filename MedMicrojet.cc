#include <vector>

#include "MedMicrojet.h"

MedMicrojet::MedMicrojet(double Rjet)
{

  _moment.clear();
  _spectrum.clear();
  _Rjet = Rjet;

}

void MedMicrojet::add_spectrum_bin(double Ejet, double moment)
{

  _spectrum.push_back( std::make_pair( Ejet, moment ) );

}
