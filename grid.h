#ifndef GRID_H
#define GRID_H

#include <vector>

typedef struct grid
{

  std::vector<double> pval, deltap;
  std::vector<double> yval, deltay;

  std::string mode;

} Grid;

#endif // Grid_H
