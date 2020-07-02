#ifndef VacMicrojet_H
#define VacMicrojet_H

#include <vector>

class VacMicrojet
{

  public:

    VacMicrojet(double Rjet);
    VacMicrojet();
    ~VacMicrojet() {}

    double moment() {return _moment;}
    void set_moment(double in_moment) {_moment = in_moment;}
    void reset_moment() {_moment = 1.0;}

    double Rjet() {return _Rjet;}

    void add_spectrum_bin(double Ejet, double spec);

    std::vector< std::pair <double,double> > spectrum() {return _spectrum;}    

  private: 

    double _moment;
    std::vector< std::pair <double,double> > _spectrum;
    double _Rjet;

};

#endif // VacMicrojet_H
