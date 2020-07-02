#ifndef cDGLAP_H
#define cDGLAP_H

#include "data.h"
#include "grid.h"
#include "ispec.h"
#include "Collimator.h"
#include "VacMicrojet.h"
#include "MedMicrojet.h"

class cDGLAP {

  private:
    
    InitData DATA;

    Grid GRID;

    InitSpec ISPEC;

    std::vector< std::pair <VacMicrojet,VacMicrojet> > vac_spectra;
    std::vector< std::pair <VacMicrojet,VacMicrojet> > med_spectra;
    std::vector< std::vector< std::pair <double,std::vector<double>> > > collim_ini;
    
    std::string mode;

  public:

    cDGLAP(std::string input_file);  
    ~cDGLAP() {}

    // Set mode
    void set_mode(std::string in_mode) { mode = in_mode; }

    // Read parameters
    InitData read_in_parameters(std::string input_file);
    InitSpec read_in_initial_spectrum();

    // Initialization
    void init();

    // Set grid
    void set_grid(double Ejet, double Rjet, std::string mode);

    // Get slope
    double slope(double Ejet, std::string flav, std::string vac_or_med);

    // Vac spectra methods
    void run_moments_vac();
    void evolve_vac_microjet(int k, VacMicrojet &VMq, VacMicrojet &VMg);
    double IntegrateVacMoment(std::string splitting, double Ejet, double theta);
    void NextVacMoment(double q_moment, double g_moment, int k, double step, double &FQ, double &FG);

    // Med spectra methods
    void run_moments_med();
    void evolve_collim(int k, Collimator &Cg, Collimator &Cq, double step_fac);
    void IntegrateCollim(std::vector<double> G_dist, std::vector<double> Q_dist, int k, double Gqw, double Qqw, double step, std::vector<double> &G, std::vector<double> &Q);
    void evolve_med_microjet(int k, MedMicrojet &MMq, MedMicrojet &MMg, Collimator &Cq, Collimator &Cg);
    void NextMedMoment(std::vector<double> MQ_dist, std::vector<double> MG_dist, int k, double step, std::vector<double> &next_MQ, std::vector<double> &next_MG, std::vector<double> yq, std::vector<double> yg);

    // Print to file
    void print_spectra();

};

#endif // cDGLAP_H
