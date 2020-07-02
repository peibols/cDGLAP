#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <iostream>

#include "Collimator.h"
#include "Util.h"

using namespace Util;

Collimator::Collimator(const InitData &DATA_in, double Ejet, double Rjet, double nspec, std::string flav, int in_nu, int in_nr) : DATA(DATA_in)
{

  if ( flav != "quark" && flav != "gluon" ) { std::cout << " Unrecognized flav in collimator constructor = " << flav << std::endl; exit(0); }
  else _flav = flav;

  if ( DATA.constant_qw ) {

    _bare_qw = DATA.bare_qw;  

  }
  else {
    
    set_bare_qw( Ejet, Rjet, nspec );  
  
  }

  _Rjet = Rjet;

  int np = DATA.np;
  for (int ip=0; ip<np+1; ip++) {
    _collim.push_back( _bare_qw );
    //_collim.push_back( 1. );
  }   
  for (int ip=0; ip<np+1; ip++) _temp_collim.push_back( 0. ); 

}

struct nloint_params
{
  double x, Q2;
  gsl_complex Omega;
  double length;
};

inline double nloint( double s, void *params )
{

  struct nloint_params *p = (struct nloint_params*) params;
  
  double x = p->x;
  double Q2 = p->Q2;
  gsl_complex Omega = p->Omega;
  double length = p->length;

  gsl_complex iw = gsl_complex_rect(0.,-0.5*x);
  gsl_complex k2;
  if (fabs(GSL_REAL(Omega))==0. && fabs(GSL_IMAG(Omega))==0. ) {
    k2 = gsl_complex_div_real( iw, -s );
  }
  else {
    gsl_complex onepiece = gsl_complex_mul( iw, Omega );
    gsl_complex twopiece = gsl_complex_sub( gsl_complex_cot( gsl_complex_mul_real(Omega,s) ) , gsl_complex_tan( gsl_complex_mul_real(Omega,length-s) ) );
    k2 = gsl_complex_mul_real( gsl_complex_mul( onepiece, twopiece ) , -1.0);
  }

  gsl_complex intpiece = gsl_complex_add_real( gsl_complex_log( gsl_complex_add_real( gsl_complex_div_real(k2,-Q2), eps1 ) ) , egamma );
  gsl_complex comp_result = gsl_complex_div( intpiece, k2 );

  return GSL_REAL(comp_result); 
}

struct radspec_params
{
  double nu_int;
  std::string flav;
  double omc;
  double Rjet, qhat, length;
  bool dla;
  double q0;
  bool running_alpha;
  int regime;
  double u2star, temp;
  double f_abarmed;
  double recoparam;
};

inline double findQ2( double w, double q, double u2 )
{
  double Q2iter = sqrt( q * w );
  for (unsigned a=0; a<10; a++) {
    if (Q2iter > u2) Q2iter = sqrt(q * w * log(Q2iter/u2));
  }
  if (Q2iter > u2) return Q2iter;
  else return u2;
}

inline double radspec( double x, void *params )
{

  struct radspec_params *p = (struct radspec_params*) params;
  
  double nu_int = p->nu_int;
  std::string flav = p->flav;
  double omc = p->omc;
  double Rjet = p->Rjet;
  double qhat = p->qhat;
  double length = p->length;
  bool dla = p->dla;
  double q0 = p->q0;
  bool running_alpha = p->running_alpha;
  int regime = p->regime;
  double u2star = p->u2star;
  double temp = p->temp;
  double const_abarmed = p->f_abarmed;
  double recoparam = p->recoparam;

  // Check cut-off
  if (x < temp + eps1) return 0;

  // Alpha
  double t = pow( x * qhat, 1./4. ); 
  double abarmed;
  if (running_alpha) {
    abarmed = alphas( t, dla, q0 ) / M_PI;
    if ( flav == "quark") abarmed *= CF;
    else abarmed *= CA;
  }
  else {
    abarmed = const_abarmed;
  }
 
  double Q2 = findQ2(x,qhat,u2star);
  if (isnan(Q2)) { std::cout << "Q2=" << Q2 << std::endl; exit(1); }
  double u = 0.5 * sqrt( qhat / x * log( Q2 / u2star ) );
  gsl_complex Omega = gsl_complex_rect(u,-u);
  
  // LO
  double clog = log( gsl_complex_abs( gsl_complex_cos( gsl_complex_mul_real(Omega,length) ) ) );
  double spec_LO = 2.0 / x * abarmed * clog;

  // NLO
  int iter = 1000;
  gsl_integration_workspace *nlo = gsl_integration_workspace_alloc (iter); 
  
  gsl_function G;
  G.function = &nloint;
  
  double rel_err = 1e-7;
 
  struct nloint_params params_nlo = { x, Q2, Omega, length };
  G.params = &params_nlo;

  double result, error;
  
  double s_lo = 0.;
  double s_hi = length;
  gsl_integration_qags (&G, s_lo, s_hi, 0, rel_err, iter,
                        nlo, &result, &error);

  if (isnan(result)) {
    std::cout << "nan result at regime= " << regime << std::endl;
    std::cout << " Real Omega= " << GSL_REAL(Omega) << " Im Omega= " << GSL_IMAG(Omega) << std::endl;
    exit(1);
  }
  double spec_NLO = -qhat * abarmed / 2.0 / x * result;
  
  gsl_integration_workspace_free(nlo);

  double angular;
  if (regime==0) {
    angular = ( 1.0 - exp( -nu_int * x * ( 1. - recoparam * Rjet * Rjet ) ) );
  }
  else if (regime==1) {
    angular = ( 1.0 - exp( -nu_int * x ) ) * exp ( -2.0 * x * x * Rjet * Rjet / qhat / length );
  }
  else {
    angular = ( 1.0 - exp( -nu_int * x ) ) * qhat * length / 4. / pow( x * Rjet ,2.);
  }

  return (spec_LO + spec_NLO) * angular;

}

void Collimator::set_bare_qw(double Ejet, double Rjet, double nspec) {

  double f_abarmed = DATA.abarmed * CF;
  if ( _flav == "gluon" ) f_abarmed *= CA / CF;

  double nu_int = nspec / Ejet;
  double oms = f_abarmed * f_abarmed * DATA.omc;
  double u2star = DATA.u2star;
  double temp = DATA.temp;
  double recoparam = DATA.recoveryparam;

  std::cout << " nu_int = " << nu_int << " wc= " << DATA.omc << " oms= " << oms << std::endl; 

  int iter = 1000;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (iter);
  
  gsl_function F;
  F.function = &radspec;
  
  double rel_err = 1e-7;
 
  // Turb regime 
  int regime=0;
  struct radspec_params params_soft = { nu_int, _flav, DATA.omc,
				  Rjet, DATA.qhat, DATA.length,
 				  DATA.dla, DATA.q0, DATA.running_alpha, regime,
 				  u2star, temp, f_abarmed, recoparam };
  F.params = &params_soft;

  double result_soft, error_soft;
  
  double w_lo = temp;
  double w_hi = oms;
  gsl_integration_qags (&F, w_lo, w_hi, 0, rel_err, iter,
                        w, &result_soft, &error_soft);
 
  
  std::cout << " result_soft= " << result_soft << " exp(-result_soft)= " << exp(-result_soft) << std::endl;
  
  // LPM regime
  regime=1;
  struct radspec_params params_hard = { nu_int, _flav, DATA.omc,
				  Rjet, DATA.qhat, DATA.length,
 				  DATA.dla, DATA.q0, DATA.running_alpha, regime,
 				  u2star, temp, f_abarmed, recoparam };
  F.params = &params_hard;

  double result_hard, error_hard;
  
  w_lo = oms;
  w_hi = DATA.omc;
  gsl_integration_qags (&F, w_lo, w_hi, 0, rel_err, iter,
                        w, &result_hard, &error_hard);
  
  std::cout << " result_hard= " << result_hard << " exp(-result_hard)= " << exp(-result_hard) << std::endl;
  
  // GLV regime
  regime=2;
  struct radspec_params params_glv = { nu_int, _flav, DATA.omc,
				  Rjet, DATA.qhat, DATA.length,
 				  DATA.dla, DATA.q0, DATA.running_alpha, regime,
 				  u2star, temp, f_abarmed, recoparam };
  F.params = &params_glv;

  double result_glv, error_glv;
  
  w_lo = DATA.omc;
  gsl_integration_qagiu (&F, w_lo, 0, rel_err, iter,
                        w, &result_glv, &error_glv);
  
  std::cout << " result_glv= " << result_glv << " exp(-result_glv)= " << exp(-result_glv) << std::endl;
 
  // Elastic term 
  double ehat = DATA.ehat;
  if (_flav=="quark") ehat *= CF/CA;
  double elastic = exp(-ehat * (1.-recoparam*Rjet*Rjet) * DATA.length * nu_int); 

  std::cout << " elastic = " << elastic << std::endl;

  // Final result 
  _bare_qw = exp( - (result_soft + result_hard + result_glv) );
  _bare_qw *= elastic;
  
  //std::cout << " result= " << result << " error= " << error << std::endl;
  //exit(0);   

  gsl_integration_workspace_free(w);
}
