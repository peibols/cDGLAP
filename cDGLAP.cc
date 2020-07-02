#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <vector>

#include "Util.h"
#include "cDGLAP.h"

cDGLAP::cDGLAP(std::string input_file) {
  DATA = read_in_parameters(input_file);
  ISPEC = read_in_initial_spectrum();
}

void cDGLAP::init() {

  vac_spectra.clear();
  med_spectra.clear();
  collim_ini.clear();

  return;
}

void cDGLAP::set_grid(double Ejet, double Rjet, std::string mode) {

  // Linear Grid for p
  
  int np = DATA.np;
  double Emin = DATA.Emin;
  
  std::vector<double> pval, deltap;
  
  for (int ip=0; ip<np+1; ip++) {
    pval.push_back( Emin * pow( Ejet / Emin , double(ip)/double(np)) );
    deltap.push_back( pval[pval.size()-1] * ( pow( Ejet / Emin , 1.0 / double(np)) - 1.0) );
  }
  
  // Log grid for y
  
  int ny = DATA.ny;
  
  std::vector<double> yval, deltay;
  double ymin, ymax;

  if ( mode == "collim_ini" ) {
  
    double thetac = DATA.thetac;
    //std::cout << " theta_c = " << thetac << std::endl;
 
    //thetac=sqrt(8./Ejet/DATA.length);
    ymin = log( Rjet / thetac );
    ymax = 0.0;
 
    for (int iy=0; iy<ny+1; iy++) { 
      deltay.push_back( - ( ymax - ymin ) / double(ny) );
      yval.push_back( ymin - double(iy) * deltay[deltay.size()-1] );
      //std::cout << "mode= " << mode << " yval= " << yval[iy] << " deltay= " << deltay[iy] << " theta= " << Rjet * exp(-yval[iy]) << std::endl;
    }
  
  }
  else if ( mode == "moments_vac" ) {

    double thetamax = DATA.thetamax;

    ymin = log( thetamax / Rjet );
    ymax = 0.0;
    
    for (int iy=0; iy<ny+1; iy++) { 
      deltay.push_back( - ( ymax - ymin ) / double(ny) );
      yval.push_back( ymin - double(iy) * deltay[deltay.size()-1] );
      //std::cout << "mode= " << mode << " yval= " << yval[iy] << " deltay= " << deltay[iy] << " theta= " << Rjet * exp(yval[iy]) << std::endl;
    }

  }
  else if ( mode == "moments_med" ) {

    double thetamax = DATA.thetamax;

    ymin = log( thetamax / Rjet );
    ymax = 0.0;

    for (int iy=0; iy<ny+1; iy++) { 
      deltay.push_back( - ( ymax - ymin ) / double(ny) );
      yval.push_back( ymin - double(iy) * deltay[deltay.size()-1] );
      //std::cout << "mode= " << mode << " yval= " << yval[iy] << " deltay= " << deltay[iy] << " theta= " << thetamax * exp(-yval[iy]) << std::endl;
    }
  
  }
  else {

    std::cout << " Unrecognized grid mode= " << mode << std::endl;
    exit(0);

  } 
   
  // Check p grid
  for (int ip=0; ip<np+1; ip++) {
    //std::cout << " pval= " << pval[ip] << " deltap= " << deltap[ip] << std::endl;
  }


  Grid new_GRID;

  new_GRID.pval = pval;    
  new_GRID.deltap = deltap;    
  new_GRID.yval = yval;    
  new_GRID.deltay = deltay;
 
  new_GRID.mode = mode;

  GRID = new_GRID;   

}

void cDGLAP::print_spectra()
{

  int maxR=DATA.binR;
  int maxP=DATA.binP;

  // Vacuum
  for (int iR=0; iR<maxR+1; iR++) {

    if ( vac_spectra.empty() ) continue;

    double thetajet = vac_spectra[iR].first.Rjet();
    std::ostringstream fvs;
    fvs << "./vac_spectra/spec_R" << thetajet*10. << ".dat";
    std::ofstream vacfile(fvs.str().c_str(),std::ios_base::binary);

    for (int iP=0; iP<maxP+1; iP++) {

      double ptbin = vac_spectra[iR].first.spectrum()[iP].first;

      double quark_cont = vac_spectra[iR].first.spectrum()[iP].second;
      double gluon_cont = vac_spectra[iR].second.spectrum()[iP].second;

      vacfile << ptbin << " " << quark_cont << " " << gluon_cont << " " << quark_cont+gluon_cont << std::endl;
    }

    vacfile.close();

  }

  // Collimator Ini 
  for (int iR=0; iR<maxR+1; iR++) {

    if ( collim_ini.empty() ) continue;

    double thetajet = vac_spectra[iR].first.Rjet();
    std::ostringstream fcs;
    fcs << "./collim_ini/collim_R" << thetajet*10. << ".dat";
    std::ofstream colfile(fcs.str().c_str(),std::ios_base::binary);
  
    for (int iP=0; iP<maxP+1; iP++) {
  
      double ptbin = collim_ini[iR][iP].first;

      double bare_q = collim_ini[iR][iP].second[0];
      double resum_q = collim_ini[iR][iP].second[1];
      double bare_g = collim_ini[iR][iP].second[2];
      double resum_g = collim_ini[iR][iP].second[3];

      colfile << ptbin << " " << bare_q << " " << resum_q << " " << bare_g << " " << resum_g << std::endl;
    }

    colfile.close();
	  
  }	  
	  
  // Medium
  for (int iR=0; iR<maxR+1; iR++) {

    if ( med_spectra.empty() ) continue;
    
    double thetajet = med_spectra[iR].first.Rjet();
    std::ostringstream fms;
    fms << "./med_spectra/spec_R" << thetajet*10. << ".dat";
    std::ofstream medfile(fms.str().c_str(),std::ios_base::binary);
  
    for (int iP=0; iP<maxP+1; iP++) {
  
      double ptbin = med_spectra[iR].first.spectrum()[iP].first;

      double quark_cont = med_spectra[iR].first.spectrum()[iP].second;
      double gluon_cont = med_spectra[iR].second.spectrum()[iP].second;

      medfile << ptbin << " " << quark_cont << " " << gluon_cont << " " << quark_cont+gluon_cont << std::endl;
    }

    medfile.close();

  }

  // RAA
  for (int iR=0; iR<maxR+1; iR++) {
    
    if ( vac_spectra.empty() || med_spectra.empty() ) continue;

    double thetajet = med_spectra[iR].first.Rjet();
    std::ostringstream frs;
    frs << "./RAA/R" << thetajet*10. << ".dat";
    std::ofstream raafile(frs.str().c_str(),std::ios_base::binary);
  
    for (int iP=0; iP<maxP+1; iP++) {
  
      double ptbin = med_spectra[iR].first.spectrum()[iP].first;

      double v_quark_cont = vac_spectra[iR].first.spectrum()[iP].second;
      double v_gluon_cont = vac_spectra[iR].second.spectrum()[iP].second;
      
      double m_quark_cont = med_spectra[iR].first.spectrum()[iP].second;
      double m_gluon_cont = med_spectra[iR].second.spectrum()[iP].second;

      raafile << ptbin << " " << m_quark_cont/v_quark_cont << " " << m_gluon_cont/v_gluon_cont << " " << (m_quark_cont+m_gluon_cont)/(v_quark_cont+v_gluon_cont) << std::endl;
    }

    raafile.close();

  }

}

double cDGLAP::slope(double Ejet, std::string flav, std::string vac_or_med)
{

  //return 5.;
  if ( flav == "quark" ) {

    if ( vac_or_med == "vac" ) return ISPEC.VQc + ISPEC.VQd * log ( Ejet / ISPEC.VQb ) + ISPEC.VQe * pow ( log ( Ejet / ISPEC.VQb ) , 2.0 ) + ISPEC.VQf * pow ( log ( Ejet / ISPEC.VQb ) , 3.0 );
    else if ( vac_or_med == "med") return ISPEC.MQc + ISPEC.MQd * log ( Ejet / ISPEC.MQb ) + ISPEC.MQe * pow ( log ( Ejet / ISPEC.MQb ) , 2.0 ) + ISPEC.MQf * pow ( log ( Ejet / ISPEC.MQb ) , 3.0 );
    else { std::cout << " Unrecognized vac_or_med in slope = " << vac_or_med << std::endl; exit(0); }

  }
  else if ( flav == "gluon" ) {
    
    if ( vac_or_med == "vac" ) return ISPEC.VGc + ISPEC.VGd * log ( Ejet / ISPEC.VGb ) + ISPEC.VGe * pow ( log ( Ejet / ISPEC.VGb ) , 2.0 ) + ISPEC.VGf * pow ( log ( Ejet / ISPEC.VGb ) , 3.0 );
    else if ( vac_or_med == "med") return ISPEC.MGc + ISPEC.MGd * log ( Ejet / ISPEC.MGb ) + ISPEC.MGe * pow ( log ( Ejet / ISPEC.MGb ) , 2.0 ) + ISPEC.MGf * pow ( log ( Ejet / ISPEC.MGb ) , 3.0 );
    else { std::cout << " Unrecognized vac_or_med in slope = " << vac_or_med << std::endl; exit(0); }

  }
  else {

    std::cout << " Unrecognized flav in slope = " << flav << std::endl;
    exit(0);

  }

}

InitData cDGLAP::read_in_parameters(std::string input_file) {

  InitData parameter_list;

  std::string tempinput;

  // Jet pt min
  double temp_pmin = 10.0;
  tempinput = Util::StringFind4(input_file, "pmin");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_pmin;
  parameter_list.pmin = temp_pmin;

  // Jet pt max
  double temp_pmax = 1000.0;
  tempinput = Util::StringFind4(input_file, "pmax");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_pmax;
  parameter_list.pmax = temp_pmax;
  
  // Switch for DLA
  bool temp_dla = false;
  tempinput = Util::StringFind4(input_file, "dla");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_dla;
  parameter_list.dla = temp_dla; 

  // Switch for DLA
  bool temp_lin_collim = false;
  tempinput = Util::StringFind4(input_file, "lin_collim");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_lin_collim;
  parameter_list.lin_collim = temp_lin_collim; 
  
  // Switch for Linear Evolution
  bool temp_linear_evol = false;
  tempinput = Util::StringFind4(input_file, "linear_evol");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_linear_evol;
  parameter_list.linear_evol = temp_linear_evol; 
  
  // Switch for Constant QW
  bool temp_constant_qw = false;
  tempinput = Util::StringFind4(input_file, "constant_qw");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_constant_qw;
  parameter_list.constant_qw = temp_constant_qw; 
  
  // Switch for Running Alpha
  bool temp_running_alpha = true;
  tempinput = Util::StringFind4(input_file, "running_alpha");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_running_alpha;
  parameter_list.running_alpha = temp_running_alpha; 
  
  // Jet pt max
  double temp_bare_qw = 0.5;
  tempinput = Util::StringFind4(input_file, "bare_qw");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_bare_qw;
  parameter_list.bare_qw = temp_bare_qw;
  
  // Temp
  double temp_temp = 0.4;
  tempinput = Util::StringFind4(input_file, "temp");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_temp;
  parameter_list.temp = temp_temp;
  
  // Temp^2
  double temp2 = 0.4*0.4;
  tempinput = Util::StringFind4(input_file, "temp2");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp2;
  
  // ehat
  double temp_ehat = 1.0/0.4;
  tempinput = Util::StringFind4(input_file, "ehat");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_ehat;
  
  // g
  double g = 1.95;
  tempinput = Util::StringFind4(input_file, "g");
  if (tempinput != "empty") std::istringstream(tempinput) >> g;
  
  // Qhat
  double temp_qhat = 1.0;
  tempinput = Util::StringFind4(input_file, "qhat");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_qhat;
  
  // QhatL^2
  double qhatL2 = 1.0;
  tempinput = Util::StringFind4(input_file, "qhatL2");
  if (tempinput != "empty") std::istringstream(tempinput) >> qhatL2;
  
  // QhatL^3
  double qhatL3 = 1.0;
  tempinput = Util::StringFind4(input_file, "qhatL3");
  if (tempinput != "empty") std::istringstream(tempinput) >> qhatL3;
  
  // qhatlogfac
  double qhatlogfac = 1.0;
  tempinput = Util::StringFind4(input_file, "qhatlogfac");
  if (tempinput != "empty") std::istringstream(tempinput) >> qhatlogfac;
  
  // Length
  double temp_length = 3.5/0.2;
  tempinput = Util::StringFind4(input_file, "length");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_length;
  parameter_list.length = temp_length;
  
  // Length
  double temp_recoveryparam = 4./M_PI/M_PI;
  tempinput = Util::StringFind4(input_file, "recoveryparam");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_recoveryparam;
  parameter_list.recoveryparam = temp_recoveryparam;
  
  // Grid size, momentum
  int temp_np = 2000;
  tempinput = Util::StringFind4(input_file, "np");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_np;
  parameter_list.np = temp_np;

  // Grid size, angle
  int temp_ny = 50;
  tempinput = Util::StringFind4(input_file, "ny");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_ny;
  parameter_list.ny = temp_ny;
 
  // LHC or RHIC
  std::string temp_system = "LHC";
  tempinput = Util::StringFind4(input_file, "system");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_system;
  parameter_list.system = temp_system;
  
  // Fill rest by hand
  parameter_list.binP = 10;
  parameter_list.binR = 8;

  parameter_list.Emin = 1.0;
  
  parameter_list.thetamin = 0.2;
  parameter_list.thetamax = 1.0;
 
  parameter_list.q0 = 0.09;
  parameter_list.lqcd = 0.15;
 
  // Medium params
  double nf_therm=3.;
  double mD = sqrt((1.+nf_therm/6.) * g * g * temp2);
  parameter_list.u2star = 0.25 * mD * mD * exp(-2.0 + 2.0*egamma);

  double ndens = (1.+nf_therm/6.);
  double as = g * g / 4. / M_PI;
  double qhat0 = 4. * M_PI * as * as * Nc * ndens;
  qhat0 *= qhatlogfac;

  double qhat = qhat0 * temp_qhat;
  parameter_list.qhat = qhat;
  std::cout << " qhat= " << qhat << std::endl;

  double ehat = temp_ehat * qhat0;
  std::cout << "ehat= " << ehat << std::endl;
  parameter_list.ehat = ehat;

  parameter_list.abarmed = as / M_PI;

  parameter_list.thetac = sqrt( 12.0 / ( qhatL3 * qhat0 ) );
  std::cout << " thetac= " <<  parameter_list.thetac << std::endl; 
  parameter_list.omc = qhatL2 * qhat0 / 2.0;

  std::cout << "u2star= " << parameter_list.u2star << std::endl;
  
  return parameter_list;

}

InitSpec cDGLAP::read_in_initial_spectrum() {

  InitSpec ispec_list;

  std::ifstream tempfile;

  // Quark Vac
  if (DATA.system == "LHC") tempfile.open("./spec_param_lhc/vac_quark.dat");
  else tempfile.open("./spec_param_rhic/vac_quark.dat");
  if (tempfile.fail()) { std::cout << " File quark_vac not found! " << std::endl; exit(0); }
  tempfile >> ispec_list.VQa >> ispec_list.VQb >> ispec_list.VQc >> ispec_list.VQd >> ispec_list.VQe >> ispec_list.VQf;
  tempfile.close();

  // Quark Med
  if (DATA.system == "LHC") tempfile.open("./spec_param_lhc/med_quark.dat");
  //if (DATA.system == "LHC") tempfile.open("./spec_param_lhc/vac_quark.dat");
  else tempfile.open("./spec_param_rhic/med_quark.dat");
  if (tempfile.fail()) { std::cout << " File quark_med not found! " << std::endl; exit(0); }
  tempfile >> ispec_list.MQa >> ispec_list.MQb >> ispec_list.MQc >> ispec_list.MQd >> ispec_list.MQe >> ispec_list.MQf;
  tempfile.close();
  
  // Gluon Vac
  if (DATA.system == "LHC") tempfile.open("./spec_param_lhc/vac_gluon.dat");
  else tempfile.open("./spec_param_rhic/vac_gluon.dat");
  if (tempfile.fail()) { std::cout << " File gluon_vac not found! " << std::endl; exit(0); }
  tempfile >> ispec_list.VGa >> ispec_list.VGb >> ispec_list.VGc >> ispec_list.VGd >> ispec_list.VGe >> ispec_list.VGf;
  tempfile.close();
  
  // Gluon Med
  if (DATA.system == "LHC") tempfile.open("./spec_param_lhc/med_gluon.dat");
  //if (DATA.system == "LHC") tempfile.open("./spec_param_lhc/vac_gluon.dat");
  else tempfile.open("./spec_param_rhic/med_gluon.dat");
  if (tempfile.fail()) { std::cout << " File gluon_med not found! " << std::endl; exit(0); }
  tempfile >> ispec_list.MGa >> ispec_list.MGb >> ispec_list.MGc >> ispec_list.MGd >> ispec_list.MGe >> ispec_list.MGf;
  tempfile.close();
  
  return ispec_list;  

}
