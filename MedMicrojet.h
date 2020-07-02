#ifndef MedMicrojet_H
#define MedMicrojet_H

#include <vector>

class MedMicrojet
{

  public:

    MedMicrojet(double Rjet);
    MedMicrojet();
    ~MedMicrojet() {}

    std::vector<double> moment() {return _moment;}
    void set_moment(std::vector<double> in_moment) {_moment = in_moment;}

    double Rjet() {return _Rjet;}

    void add_spectrum_bin(double Ejet, double spec);

    std::vector< std::pair <double,double> > spectrum() {return _spectrum;}    

  private: 

    std::vector<double> _moment;
    std::vector< std::pair <double,double> > _spectrum;
    double _Rjet;

};

#endif // MedMicrojet_H
