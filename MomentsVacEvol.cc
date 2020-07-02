#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "Util.h"
#include "cDGLAP.h"
#include "VacMicrojet.h"

using namespace Util;

void cDGLAP::run_moments_vac() {

  set_mode("vac");

  double pmin=DATA.pmin;
  double pmax=DATA.pmax;

  double thetamin=DATA.thetamin;
  double thetamax=DATA.thetamax;

  int maxP=DATA.binP;
  int maxR=DATA.binR;

  // Grid size
  int ny = DATA.ny;

  auto pp = new double[maxP+1]; 
  for (int iP=0; iP<maxP+1; iP++) {
    double pjet = pmin * pow( pmax / pmin, double(iP) / double(maxP) );
    pp[iP] = pjet;
  }

  // Evolution from thetamax to thetajet
  for (int iR=0; iR<maxR+1; iR++) {

    double thetajet = thetamin + (thetamax - thetamin) * double (iR) / double(maxR);

    VacMicrojet VMq( thetajet );
    VacMicrojet VMg( thetajet );
    
    for (int iP=0; iP<maxP+1; iP++) {
      
      // Check range of parametrisation
      if ( pp[iP] < ISPEC.VQb || pp[iP] < ISPEC.MQb || pp[iP] < ISPEC.VGb || pp[iP] < ISPEC.MGb ) {
        std::cout << " Ejet = " << pp[iP] << " too small, cannot use this initial spectrum parametrisation " << std::endl;
        continue;
      }
      
      set_grid( pp[iP], thetajet, "moments_vac" );
      
      VMq.reset_moment();
      VMg.reset_moment();
      
      std::cout << "Jet Pt = " << pp[iP] << " Jet R = " << thetajet << std::endl;

      // Evolution in y
      for (int iy=0; iy<ny; iy++) {
        
        evolve_vac_microjet( iy, VMq, VMg );

      } // End y loop
      
      std::cout << " Final Vac Q moment = " << VMq.moment() << std::endl;
      std::cout << " Final Vac G moment = " << VMg.moment() << std::endl;

      double quark_ispec = ISPEC.VQa * pow( ISPEC.VQb / pp[iP] , slope( pp[iP], "quark", "vac" ) );
      //quark_ispec=1.;
      VMq.add_spectrum_bin( pp[iP], quark_ispec * VMq.moment() );
      
      double gluon_ispec = ISPEC.VGa * pow( ISPEC.VGb / pp[iP] , slope( pp[iP], "gluon", "vac" ) );
      //gluon_ispec=1.;
      VMg.add_spectrum_bin( pp[iP], gluon_ispec * VMg.moment() );

    } // End jet momentum loop

    // Store for future use
    vac_spectra.push_back( std::make_pair( VMq, VMg ) );

  } // End jet radius loop

}

void cDGLAP::evolve_vac_microjet(int k, VacMicrojet &VMq, VacMicrojet &VMg)
{

  double q_moment = VMq.moment();
  double g_moment = VMg.moment();

  std::vector<double> deltay = GRID.deltay;

  // Linear

  if ( DATA.linear_evol ) {
    double step = deltay[k];
    double Q1=0.;
    double G1=0.;
    NextVacMoment( q_moment, g_moment, k, step, Q1, G1 );
    double q_moment_next_step = q_moment + step * Q1;
    double g_moment_next_step = g_moment + step * G1;
    VMq.set_moment(q_moment_next_step);
    VMg.set_moment(g_moment_next_step);
    return;
  }
  
  // RK4

  // First integral
  double step = 0.; 
  double Q1 = 0.;
  double G1 = 0.;
  NextVacMoment( q_moment, g_moment, k, step, Q1, G1 );

  // Second integral
  step = deltay[k] / 2.;
  double iQ = q_moment + step * Q1;
  double iG = g_moment + step * G1;
  double Q2 = 0.;
  double G2 = 0.;
  NextVacMoment( iQ, iG, k, step, Q2, G2 );

  // Third integral
  step = deltay[k] / 2.;
  iQ = q_moment + step * Q2;
  iG = g_moment + step * G2;
  double Q3 = 0.;
  double G3 = 0.;
  NextVacMoment( iQ, iG, k, step, Q3, G3 );

  // Fourth integral
  step = deltay[k];
  iQ = q_moment + step * Q3;
  iG = g_moment + step * G3;
  double Q4 = 0.;
  double G4 = 0.;
  NextVacMoment( iQ, iG, k, step, Q4, G4 );

  // Weighted sum

  double q_moment_next_step = q_moment + deltay[k] / 6.0 * 
	  		( Q1 + 2.0 * ( Q2 + Q3 ) + Q4 );
  double g_moment_next_step = g_moment + deltay[k] / 6.0 * 
	  		( G1 + 2.0 * ( G2 + G3 ) + G4 );
  VMq.set_moment(q_moment_next_step);
  VMg.set_moment(g_moment_next_step);
  
  return;

}

void cDGLAP::NextVacMoment(double q_moment, double g_moment, int k, double step, double &FQ, double &FG)
{

  std::vector<double> yval = GRID.yval;
  std::vector<double> pval = GRID.pval; // Do I need?

  double Ejet = pval[pval.size()-1];   // Do I need?
  double Rjet = DATA.thetamax * exp( -yval[0] );

  double theta = Rjet * exp( yval[k] - step );
  //std::cout << " Ejet= " << Ejet << " Rjet= " << Rjet << " theta= " << theta << std::endl;

  double GammaQQ = IntegrateVacMoment( "qq", Ejet, theta );
  double GammaGG = IntegrateVacMoment( "gg", Ejet, theta );
  double GammaGQ = IntegrateVacMoment( "gq", Ejet, theta );
  double GammaQG = IntegrateVacMoment( "qg", Ejet, theta );
  
  FQ = GammaQQ * q_moment + GammaGQ * g_moment;
  FG = GammaGG * g_moment + GammaQG * q_moment;

  return; 

}

struct mom_vac_params
{
  double nspec;
  double Ejet, theta;
  std::string splitting;
  bool dla;
  double q0;
};

inline double mom_vac( double x, void *params)
{

  struct mom_vac_params *p = (struct mom_vac_params*) params;

  double nspec = p->nspec;
  double Ejet = p->Ejet;
  double theta= p->theta;
  std::string splitting = p->splitting;
  bool dla = p->dla;
  double q0 = p->q0;

  double kt = x * ( 1.0 - x ) * Ejet * theta;

  double integrand;

  if ( splitting == "qq" ) {
    integrand = Pqq( x, dla ) * ( pow( x, nspec - 1.0 ) - 1.0 );
  }
  else if ( splitting == "gg" ) {
    integrand = Pgg( x, dla ) * ( pow( x, nspec - 1.0 ) - x ) - Pqg( x, dla );
  }
  else if ( splitting == "gq" ) {
    integrand = Pgq( x, dla ) * pow( x, nspec - 1.0 );
  }
  else if ( splitting == "qg" ) {
    integrand = Pqg( x, dla ) * pow( x, nspec - 1.0 );
  }
  else {
    std::cout << " Unrecognized splitting kernel = " << splitting << std::endl;
    exit(0);
  }

  //std::cout << " integrand= " << integrand << " alphas= " << alphas( kt, dla, q0 ) << std::endl;
  return integrand * alphas( kt, dla, q0 ) / M_PI;

}

double cDGLAP::IntegrateVacMoment(std::string splitting, double Ejet, double theta)
{

  std::string flav = "quark";
  if ( splitting == "gg" || splitting == "qg" ) flav = "gluon";
  double nspec = slope( Ejet, flav, mode );

  int iter = 1000;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (iter);
 
  gsl_set_error_handler_off();
 
  gsl_function F;
  F.function = &mom_vac;
  
  double rel_err = 1e-7;
  
  bool dla = DATA.dla;
  double q0 = DATA.q0;

  struct mom_vac_params params = { nspec, Ejet, theta, splitting, dla, q0 };
  F.params = &params;

  double result, error;
  
  double z_lo = 0.0;
  double z_hi = 1.0;

  //std::cout << " Doing splitting = " << splitting << " Ejet= " << Ejet << " theta= " << theta << std::endl;
  int status; 
  do {
   
    status = gsl_integration_qags (&F, z_lo, z_hi, 0, rel_err, iter,
                        w, &result, &error);
    
    rel_err *= 1.1;
  
  } while ( status == GSL_EROUND );

  if (rel_err> 1e-6 ) std::cout << " result= " << result << " error= " << error << std::endl;

  gsl_integration_workspace_free(w);

  return result;
}
