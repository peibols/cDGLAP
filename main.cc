#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <assert.h>
#include <fstream>

#include "cDGLAP.h"

int main(int argc, char **argv)
{

  assert(argc==2);
 
  std::string input_file;
  input_file = *(argv+1);

  std::cout << "Start Program" << std::endl;

  cDGLAP cdglap(input_file); 

  cdglap.init();

  // Vacuum spectra
  cdglap.run_moments_vac();

  // Medium spectra
  cdglap.run_moments_med();

  // Print spectra and RAA to file
  cdglap.print_spectra();

  return 0;

  std::cout << "End Program" << std::endl;

}
