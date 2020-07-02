#ifndef DATA_H
#define DATA_H

typedef struct init_data {

  // Range of ptjet
  double pmin, pmax;
  int binP;
  double Emin;

  // Range of Rjet
  double thetamin, thetamax;
  int binR;   

  // Grid size
  int np;
  int ny;

  // Physical parameters
  double q0;
  double lqcd;

  // Medium params
  double qhat;
  double ehat;
  double length;
  double temp;
  double u2star;
  double abarmed;
  double thetac;
  double omc;
  double bare_qw;		// Used if considered constant

  // Uncertainties
  double recoveryparam;

  // Switches
  bool dla;
  bool linear_evol;
  bool constant_qw;
  bool lin_collim;
  bool running_alpha;

  // System
  std::string system;

} InitData;

#endif // DATA_H
