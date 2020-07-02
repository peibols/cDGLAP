#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <omp.h>

#include <iostream>
#include <cmath>

#include "cDGLAP.h"
#include "Util.h"

using namespace std;
using namespace Util;

void cDGLAP::run_moments_med() {

  set_mode("med");
  
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

  for (int iR=0; iR<maxR+1; iR++) {

    bool unresolved_jet=0;
  
    double thetajet = thetamin + (thetamax - thetamin) * double (iR) / double(maxR);

    VacMicrojet MMq( thetajet );
    VacMicrojet MMg( thetajet );

    std::vector< std::pair <double,std::vector<double>> > temp_collim_ini;

    for (int iP=0; iP<maxP+1; iP++) {
      
      // Check range of parametrisation
      if ( pp[iP] < ISPEC.VQb || pp[iP] < ISPEC.MQb || pp[iP] < ISPEC.VGb || pp[iP] < ISPEC.MGb ) {
        std::cout << " Ejet = " << pp[iP] << " too small, cannot use this initial spectrum parametrisation " << std::endl;
        continue;
      }

      double thisthetac=DATA.thetac;
      //thisthetac=sqrt(8./pp[iP]/DATA.length);
      if (thetajet < thisthetac) {
        std::cout << "Unresolved Rjet = " << thetajet << ", smaller than thetac = " << DATA.thetac << std::endl;
        unresolved_jet=1;
        //continue;
      }
      std::cout << "Jet Pt = " << pp[iP] << " Jet R = " << thetajet << std::endl;

      std::cout << " Collimator Ini \n";
      // Collimator initial condition
      if (!unresolved_jet) set_grid( pp[iP], thetajet, "collim_ini" );

      int in_nr = iR;
      int in_nu = iP;

      double nspec_g = slope( pp[iP], "gluon", "med" );
      Collimator Cg( DATA , pp[iP], thetajet, nspec_g, "gluon", in_nu, in_nr );

      double nspec_q = slope( pp[iP], "quark", "med" );
      Collimator Cq( DATA , pp[iP], thetajet, nspec_q, "quark", in_nu, in_nr );
     

       if (!unresolved_jet) { 
        // Evolution in y for collimator initial condition
        for (int iy=0; iy<ny; iy++) {
        
          evolve_collim( iy, Cg, Cq, 1.0 );
        
        } // End y loop
      }

      // Store Collim Ini
      vector<double> collims;
      collims.push_back( Cq.bare_qw() );
      collims.push_back( Cq.collim()[Cq.collim().size()-1] );
      collims.push_back( Cg.bare_qw() );
      collims.push_back( Cg.collim()[Cg.collim().size()-1] );
      temp_collim_ini.push_back( std::make_pair( pp[iP], collims ) );

      std::cout << " bare Q= " << collims[0] << " resum Q= " << collims[1] << " bare G= " << collims[2] << " resum G= " << collims[3] << std::endl;

      //exit(0);
      //continue;

      std::cout << " Moment Med \n";
      // Med Moment Evolution
      set_grid( pp[iP], thetajet, "moments_med" );
     
      // Initial condition 
      MMq.set_moment( Cq.collim()[Cq.collim().size()-1] );      
      MMg.set_moment( Cg.collim()[Cg.collim().size()-1] );
      //MMq.set_moment( 1.0 );      
      //MMg.set_moment( 1.0 );

      // Evolve microjet in medium
      for (int iy=0; iy<ny; iy++) {

	evolve_vac_microjet( iy, MMq, MMg );
        //evolve_med_microjet( iy, MMq, MMg, Cq, Cg );

      }      
      
      std::cout << " Final Med Q moment = " << MMq.moment() << std::endl;
      std::cout << " Final Med G moment = " << MMg.moment() << std::endl;

      double quark_ispec = ISPEC.MQa * pow( ISPEC.MQb / pp[iP] , slope( pp[iP], "quark", "med" ) );
      //quark_ispec=1.;
      MMq.add_spectrum_bin( pp[iP], quark_ispec * MMq.moment() );

      double gluon_ispec = ISPEC.MGa * pow( ISPEC.MGb / pp[iP] , slope( pp[iP], "gluon", "med" ) );
      //gluon_ispec=1.;
      MMg.add_spectrum_bin( pp[iP], gluon_ispec * MMg.moment() );

    } // End jet momentum loop

    //exit(0);
    
    // Store for future use
    med_spectra.push_back( std::make_pair( MMq, MMg ) );
    collim_ini.push_back( temp_collim_ini );

  } // End jet radius loop 

}

void cDGLAP::evolve_med_microjet(int k, MedMicrojet &MMq, MedMicrojet &MMg, Collimator &Cq, Collimator &Cg)
{

  vector<double> mq = MMq.moment();
  vector<double> mg = MMg.moment();

  std::vector<double> deltay = GRID.deltay;

  // Linear

  if ( DATA.linear_evol ) {

    // Evolve collim one step
    evolve_collim( k, Cq, Cg, 1.0 );
    vector<double> yq = Cq.collim();
    vector<double> yg = Cg.collim();

    double step = deltay[k];
    vector<double> MG1, MQ1;
    NextMedMoment( mq, mg, k, step, MQ1, MG1, yq, yg ); 
    
    vector<double> MMG_lin_next_step = mg;
    vector<double> MMQ_lin_next_step = mq;
    for (unsigned int i=0; i<mq.size(); i++) {
      MMG_lin_next_step[i] += deltay[k] * MG1[i];
      MMQ_lin_next_step[i] += deltay[k] * MQ1[i];
    }
    MMg.set_moment(MMG_lin_next_step);
    MMq.set_moment(MMQ_lin_next_step);

    return;
  }

  // RK4

  // Collimators at h = 0
  vector<double> yg_0 = Cg.collim();
  vector<double> yq_0 = Cq.collim();
  
  // First integral
  double step = 0;
  vector<double> MG1, MQ1;
  NextMedMoment( mq, mg, k, step, MQ1, MG1, yq_0, yg_0 ); 

  // Collimators at h = step / 2
  evolve_collim( k, Cq, Cg, 0.5 );
  vector<double> yq_0p5 = Cq.temp_collim();
  vector<double> yg_0p5 = Cg.temp_collim();

  // Second integral
  step = deltay[k] / 2.;
  vector<double> MG2, MQ2;
  vector<double> imq = mq;
  vector<double> img = mg;
  for (unsigned int i=0; i<mg.size(); i++) img[i] += step * MG1[i], imq[i] += step * MQ1[i];
  NextMedMoment( imq, img, k, step, MQ2, MG2, yq_0p5, yg_0p5 ); 
  
  // Third integral
  step = deltay[k] / 2.;
  vector<double> MG3, MQ3;
  imq = mq;
  img = mg;
  for (unsigned int i=0; i<mg.size(); i++) img[i] += step * MG2[i], imq[i] += step * MQ2[i];
  NextMedMoment( imq, img, k, step, MQ3, MG3, yq_0p5, yg_0p5 ); 

  // Collimators at h = step
  evolve_collim( k, Cq, Cg, 1.0 );
  vector<double> yq_1 = Cq.collim();
  vector<double> yg_1 = Cg.collim();
  
  // Fourth integral
  step = deltay[k];
  vector<double> MG4, MQ4;
  imq = mq;
  img = mg;
  for (unsigned int i=0; i<mg.size(); i++) img[i] += step * MG3[i], imq[i] += step * MQ3[i];
  NextMedMoment( imq, img, k, step, MQ4, MG4, yq_1, yg_1 ); 
  
  // Weighted sum

  vector<double> MMG_lin_next_step = mg;
  vector<double> MMQ_lin_next_step = mq;
  for (unsigned int i=0; i<mg.size(); i++) {
    MMQ_lin_next_step[i] += deltay[k] / 6.0 *
                        ( MQ1[i] + 2.0 * ( MQ2[i] + MQ3[i] ) + MQ4[i] );
    MMG_lin_next_step[i] += deltay[k] / 6.0 *
                        ( MG1[i] + 2.0 * ( MG2[i] + MG3[i] ) + MG4[i] );
  }
  MMg.set_moment(MMG_lin_next_step);
  MMq.set_moment(MMQ_lin_next_step);

  return;

}

void cDGLAP::evolve_collim(int k, Collimator &Cg, Collimator &Cq, double step_fac)
{

  vector<double> yg = Cg.collim();
  vector<double> yq = Cq.collim();

  vector<double> deltay = GRID.deltay;
  
  double Gqw = Cg.bare_qw();
  double Qqw = Cq.bare_qw();
//  cout << " doing angle k= " << k << endl;

  // Linear

  if ( DATA.linear_evol ) {
    double step = deltay[k];
    vector<double> G1, Q1;
    IntegrateCollim( yg, yq, k, Gqw, Qqw, step, G1, Q1 );
    vector<double> G_lin_next_step = yg;
    vector<double> Q_lin_next_step = yq;
    for (unsigned int i=0; i<yq.size(); i++) {
      G_lin_next_step[i] += deltay[k] * G1[i];
      Q_lin_next_step[i] += deltay[k] * Q1[i];
    }
    Cg.set_collim(G_lin_next_step);
    Cq.set_collim(Q_lin_next_step);
    return;
  }
  
  // RK4

  // First integral
  double step = 0.;
  vector<double> G1, Q1;
  IntegrateCollim( yg, yq, k, Gqw, Qqw, step, G1, Q1 );
  
  // Second integral
  step = deltay[k] / 2. * step_fac;
  vector<double> iG = yg;
  vector<double> iQ = yq;
  for (unsigned int i=0; i<yg.size(); i++) iG[i] += step * G1[i], iQ[i] += step * Q1[i];
  vector<double> G2, Q2;
  IntegrateCollim( iG, iQ, k, Gqw, Qqw, step, G2, Q2 );
  

  // Third integral
  step = deltay[k] / 2. * step_fac;
  iG = yg;
  iQ = yq;
  for (unsigned int i=0; i<yg.size(); i++) iG[i] += step * G2[i], iQ[i] += step * Q2[i];
  vector<double> G3, Q3;
  IntegrateCollim( iG, iQ, k, Gqw, Qqw, step, G3, Q3 );

  // Fourth integral
  step = deltay[k] * step_fac;
  iG = yg;
  iQ = yq;
  for (unsigned int i=0; i<yg.size(); i++) iG[i] += step * G3[i], iQ[i] += step * Q3[i];
  vector<double> G4, Q4;
  IntegrateCollim( iG, iQ, k, Gqw, Qqw, step, G4, Q4 );
  
  // Weighted sum

  vector<double> Q_next_step = yq;
  vector<double> G_next_step = yg;
  for (unsigned int i=0; i<yg.size(); i++) {
    Q_next_step[i] += deltay[k] / 6.0 *
                        ( Q1[i] + 2.0 * ( Q2[i] + Q3[i] ) + Q4[i] );
    G_next_step[i] += deltay[k] / 6.0 *
                        ( G1[i] + 2.0 * ( G2[i] + G3[i] ) + G4[i] );
  }

  if ( step_fac == 1.0 ) {
    Cg.set_collim(G_next_step);
    Cq.set_collim(Q_next_step);
  }
  else {
    Cg.set_temp_collim(G_next_step);
    Cq.set_temp_collim(Q_next_step);
  }

  return;

}

void cDGLAP::IntegrateCollim(vector<double> G_dist, vector<double> Q_dist, int k, double Gqw, double Qqw, double step, vector<double> &G, vector<double> &Q)
{

  double qhat = DATA.qhat;
  double Emin = DATA.Emin;

  int np = DATA.np;
  vector<double> pval = GRID.pval;
  vector<double> yval = GRID.yval;

  double Ejet = pval[pval.size()-1];
  
  double Rjet, theta; 

  if ( GRID.mode == "collim_ini" ) {
   
    Rjet = DATA.thetac * exp( yval[0] );
    theta = Rjet * exp( -yval[k] + step );
  
  }
  else if ( GRID.mode == "moments_med" ) {
  
    Rjet = DATA.thetamax * exp( -yval[0] );
    theta = DATA.thetamax * exp( -yval[k] + step );
  
  }
  else {

    std::cout << " Unrecognized mode for collim grid = " << GRID.mode << std::endl;
    exit(0);

  } 

  //cout << " theta= " << theta << " Rjet= " << Rjet << endl;

  double LLfac = 1.;
  double ak = LLfac * pow( 2./3. * qhat / pow( theta , 4. ), 1./3. );
  //double ak = 2./theta/theta/DATA.length;
  //ak = pow(2./3. * qhat / pow(Rjet,4.), 1./3.);

  #pragma omp parallel
  {
  
  // Output num threads
  //std::cout << " num_threads= " << omp_get_num_threads() << std::endl;

  vector<double> G_private, Q_private;

  #pragma omp for nowait schedule(static)
  for (int ip=0; ip<np+1; ip++) {
    
    double G_result=0.;
    double Q_result=0.;
    
    double w = pval[ip];

    // Phase space constraint
    if ( w < 4.0 * ak ) {
      G_private.push_back( G_result );
      Q_private.push_back( Q_result );
      continue;
    }

    double wp = w / 2.0 * ( 1.0 + sqrt(1.0 - 4.0 * ak / w) );
    double wm = w / 2.0 * ( 1.0 - sqrt(1.0 - 4.0 * ak / w) );

    wm = max ( wm , Emin );
    wp = max ( wp , Emin );
    if (wm == wp ) {
      G_private.push_back( G_result );
      Q_private.push_back( Q_result );
      continue;
    }
   
    // Construct the arrays to integrate
    double G_coll_int[np+1];
    double Q_coll_int[np+1];

    for (int jp=0; jp<np+1; jp++) {
      
      double kernel_gg = 0.;
      double kernel_qg = 0.;
      double kernel_gq = 0.;      
 
      double edif = w - pval[jp];      

      if ( edif >= Emin ) {
        
        int id = int (np * log( edif / Emin ) / log ( Ejet / Emin ));

        double z = pval[jp] / w;
        double kt = z * ( 1.0 - z ) * w * theta;

        kernel_gg = alphas( kt, DATA.dla, DATA.q0 ) / M_PI * Pgg( pval[jp] / w, DATA.dla ) / w;
	if ( DATA.lin_collim ) kernel_gg *= G_dist[ip] * ( Gqw - 1.0 );
	//else kernel_gg *= ( G_dist[jp] * G_dist[id] - G_dist[ip] );
	else kernel_gg *= ( G_dist[jp] * G_dist[id] - G_dist[ip] );
     
        kernel_qg = alphas( kt, DATA.dla, DATA.q0 ) / M_PI * Pqg( pval[jp] / w, DATA.dla ) / w;
	if ( DATA.lin_collim ) kernel_qg *= 0.;
	else kernel_qg *= ( Q_dist[jp] * Q_dist[id] - G_dist[ip] );

        kernel_gq = alphas( kt, DATA.dla, DATA.q0 ) / M_PI * Pgq( pval[jp] / w, DATA.dla ) / w;
	if ( DATA.lin_collim ) kernel_gq *= Q_dist[ip] * ( Gqw - 1.0 );
	else kernel_gq *= ( G_dist[jp] * Q_dist[id] - Q_dist[ip] );
 
      }

      G_coll_int[jp] = kernel_gg + kernel_qg;
      Q_coll_int[jp] = kernel_gq;
    
    }

    // Arrange in arrays first (needed for GSL)
    double pval_arr[pval.size()];
    copy(pval.begin(), pval.end(), pval_arr);

    // Manual integration
    /*
    double man_result=0.;
    for (unsigned int i=0; i<pval.size()-1; i++) {
      if (pval_arr[i]>=wm && pval_arr[i+1]<wp) {
        double binwidth = pval_arr[i+1]-pval_arr[i];
        double bincen = (pval_arr[i+1]+pval_arr[i])/2.;
        double binheight=(coll_int[i]+coll_int[i+1])/2.;
        double piece = binheight * binwidth;
        man_result += piece;
      }
    }
    //cout << " man_result= " << man_result << endl;
    */


    // G int    

    gsl_interp_accel *g_acc = gsl_interp_accel_alloc ();
    gsl_spline *g_spline = gsl_spline_alloc (gsl_interp_linear, np+1); 
    gsl_spline_init (g_spline, pval_arr , G_coll_int , np+1);

    G_result = gsl_spline_eval_integ(g_spline, wm, wp, g_acc);
    gsl_spline_free (g_spline);
    gsl_interp_accel_free (g_acc);

    G_private.push_back( G_result );
    
    // Q int    

    gsl_interp_accel *q_acc = gsl_interp_accel_alloc ();
    gsl_spline *q_spline = gsl_spline_alloc (gsl_interp_linear, np+1); 
    gsl_spline_init (q_spline, pval_arr , Q_coll_int , np+1);

    Q_result = gsl_spline_eval_integ(q_spline, wm, wp, q_acc);
    gsl_spline_free (q_spline);
    gsl_interp_accel_free (q_acc);

    Q_private.push_back( Q_result );

  }

  #pragma omp for schedule(static) ordered
  for (int i=0; i<omp_get_num_threads(); i++) {
    #pragma omp ordered
    G.insert(G.end(), G_private.begin(), G_private.end());
    Q.insert(Q.end(), Q_private.begin(), Q_private.end());
  }

  } //End omp parallel

}
    
//NextMedMoment( mq, mg, k, step, MQ1, MG1, yq, yg ); 
void cDGLAP::NextMedMoment(vector<double> MQ_dist, vector<double> MG_dist, int k, double step, vector<double> &next_MQ, vector<double> &next_MG, vector<double> yq, vector<double> yg)
{

  double qhat = DATA.qhat;
  double Emin = DATA.Emin;

  int np = DATA.np;
  vector<double> pval = GRID.pval;
  vector<double> yval = GRID.yval;

  double Ejet = pval[pval.size()-1];
  
  double Rjet, theta; 
  Rjet = DATA.thetamax * exp( -yval[0] );
  theta = DATA.thetamax * exp( -yval[k] + step );
  
  //cout << " theta= " << theta << " Rjet= " << Rjet << endl;

  double ak = pow( 2./3. * qhat / pow( theta , 4. ), 1./3. );
  //ak = pow(2./3. * qhat / pow(Rjet,4.), 1./3.);

  #pragma omp parallel
  {

  vector<double> MG_private, MQ_private;

  #pragma omp for nowait schedule(static)
  for (int ip=0; ip<np+1; ip++) {
    
    double MG_result=0.;
    double MQ_result=0.;
    
    double w = pval[ip];

    // Phase space constraint
    if ( w < 4.0 * ak ) {
      MG_private.push_back( MG_result );
      MQ_private.push_back( MQ_result );
      continue;
    }

    double wp = w / 2.0 * ( 1.0 + sqrt(1.0 - 4.0 * ak / w) );
    double wm = w / 2.0 * ( 1.0 - sqrt(1.0 - 4.0 * ak / w) );

    wm = max ( wm , Emin );
    wp = max ( wp , Emin );
    if (wm == wp ) {
      MG_private.push_back( MG_result );
      MQ_private.push_back( MQ_result );
      continue;
    }
   
    // Construct the arrays to integrate
    double MG_coll_int[np+1];
    double MQ_coll_int[np+1];

    for (int jp=0; jp<np+1; jp++) {
      
      double kernel_qq = 0.;
      double kernel_gg = 0.;
      double kernel_qg = 0.;
      double kernel_gq = 0.;      
 
      double edif = w - pval[jp];      

      if ( edif >= Emin ) {
        
        int id = int (np * log( edif / Emin ) / log ( Ejet / Emin ));

        double z = pval[jp] / w;
        double kt = z * ( 1.0 - z ) * w * theta;

	// Compute slopes (ignore z dependence for now)
        double nspec_q = slope( Ejet, "quark", "med" ) - 1.0;
        double nspec_g = slope( Ejet, "gluon", "med" ) - 1.0;

        kernel_gg = alphas( kt, DATA.dla, DATA.q0 ) / M_PI * Pgg( pval[jp] / w, DATA.dla ) / w;
	kernel_gg *= ( pow( z, nspec_g ) * yq[id] * MG_dist[jp] - z * MG_dist[ip] );
//	kernel_gg *= G_dist[ip] * ( Gqw - 1.0 );
     
        kernel_qg = alphas( kt, DATA.dla, DATA.q0 ) / M_PI * Pqg( pval[jp] / w, DATA.dla ) / w;
	kernel_qg *= ( pow( z, nspec_g ) * yg[id] * MQ_dist[ip] - MG_dist[ip] );

        kernel_gq = alphas( kt, DATA.dla, DATA.q0 ) / M_PI * Pgq( pval[jp] / w, DATA.dla ) / w;
	kernel_gq *= pow( z, nspec_q ) * yq[id] * MG_dist[ip];
 
        kernel_qq = alphas( kt, DATA.dla, DATA.q0 ) / M_PI * Pqq( pval[jp] / w, DATA.dla ) / w;
	kernel_qq *= ( pow( z, nspec_q ) * yg[id] * MQ_dist[jp] - MQ_dist[ip] );

      }


      MG_coll_int[jp] = kernel_gg + kernel_qg;
      MQ_coll_int[jp] = kernel_gq + kernel_qq;
    
    }

    // Arrange in arrays first (needed for GSL)
    double pval_arr[pval.size()];
    copy(pval.begin(), pval.end(), pval_arr);

    // Manual integration
    /*
    double man_result=0.;
    for (unsigned int i=0; i<pval.size()-1; i++) {
      if (pval_arr[i]>=wm && pval_arr[i+1]<wp) {
        double binwidth = pval_arr[i+1]-pval_arr[i];
        double bincen = (pval_arr[i+1]+pval_arr[i])/2.;
        double binheight=(coll_int[i]+coll_int[i+1])/2.;
        double piece = binheight * binwidth;
        man_result += piece;
      }
    }
    //cout << " man_result= " << man_result << endl;
    */

    // G int    

    gsl_interp_accel *g_acc = gsl_interp_accel_alloc ();
    gsl_spline *g_spline = gsl_spline_alloc (gsl_interp_linear, np+1); 
    gsl_spline_init (g_spline, pval_arr , MG_coll_int , np+1);

    MG_result = gsl_spline_eval_integ(g_spline, wm, wp, g_acc);
    gsl_spline_free (g_spline);
    gsl_interp_accel_free (g_acc);

    MG_private.push_back( MG_result );
    
    // Q int    

    gsl_interp_accel *q_acc = gsl_interp_accel_alloc ();
    gsl_spline *q_spline = gsl_spline_alloc (gsl_interp_linear, np+1); 
    gsl_spline_init (q_spline, pval_arr , MQ_coll_int , np+1);

    MQ_result = gsl_spline_eval_integ(q_spline, wm, wp, q_acc);
    gsl_spline_free (q_spline);
    gsl_interp_accel_free (q_acc);

    MQ_private.push_back( MQ_result );

  }
  
  #pragma omp for schedule(static) ordered
  for (int i=0; i<omp_get_num_threads(); i++) {
    #pragma omp ordered
    next_MG.insert(next_MG.end(), MG_private.begin(), MG_private.end());
    next_MQ.insert(next_MQ.end(), MQ_private.begin(), MQ_private.end());
  }

  } //End omp parallel

}
