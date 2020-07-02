#ifndef Collimator_H
#define Collimator_H

#include <vector>

#include "data.h"

class Collimator
{

  public:

    Collimator(const InitData &DATA_in, double Ejet, double Rjet, double nspec, std::string flav, int in_nu, int in_nr);
    Collimator();
    ~Collimator() {}

    std::vector<double> collim() {return _collim;}
    void set_collim(std::vector<double> in_collim) {_collim = in_collim;}

    std::vector<double> temp_collim() {return _temp_collim;}
    void set_temp_collim(std::vector<double> in_collim) {_temp_collim = in_collim;}
    
    void set_bare_qw(double Ejet, double Rjet, double nspec);
    double bare_qw() {return _bare_qw;}

  private: 

    std::vector<double> _collim;
    std::vector<double> _temp_collim;
    double _bare_qw;

    double _Rjet;

    std::string _flav;

    const InitData &DATA;
};

#endif // Collimator_H
