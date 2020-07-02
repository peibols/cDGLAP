#include "Random.h"
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <iostream>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <sstream>

#include "global.h"
#include "vector_operators.h"

#define DO_WAKE
#define JET_TOOLS

char Sfile[100];	//Source File name: one file per event

using std::vector;
using namespace std;

void read_nuclear(int, std::string);
void read_hydro(int, std::string);

void gxy(double &, double &, double &, numrand &);
double gVx(double tau, double x, double y, double eta);
double gVy(double tau, double x, double y, double eta);
double gVz(double tau, double x, double y, double eta);
double gT(double tau, double x, double y, double eta);

// Nev, Expansion, yes(1) or no(0)
int main(int argc, char** argv)
{
	assert(argc==4);
	int Nev=atoi(argv[1]);
	bool no_expansion;
	if (atoi(argv[2])==1) no_expansion=0;
	else no_expansion=1;
	std::string cent=argv[3];

	int N=1;

	std::ostringstream fq;
        fq << "./qhat_dists/cent_" << cent.c_str() << ".dat";
	std::ofstream qhatfile(fq.str().c_str(),std::ios_base::binary);
	qhatfile << "#<qhat> <qhatL^2> <qhatL^3> <L> <T> <T^2> <ehat>" << endl;	

	std::ostringstream avefq;
        avefq << "./qhat_dists/ave_cent_" << cent.c_str() << ".dat";
	std::ofstream aveqhatfile(avefq.str().c_str(),std::ios_base::binary);
	aveqhatfile << "#<qhat> <qhatL^2> <qhatL^3> <L> <T> <T^2> <ehat>" << endl;	

        double avevals[7]={0.};
	
        double step=0.02; //fm
	double Tc=0.145;
	double hbarc=0.197327;

	//Initialize Random Seed
	int njob=1;
	numrand nr(1346+njob);
	//cout << " rando= " << nr.rando() << endl;

	//Total Event Loop
	int totcount=0;
	do {

		read_nuclear(totcount+1, cent);
		read_hydro(totcount+1, cent);

		//Hydro Event Loop
		int count=0;
		double tempcent=gT(0.6,0.,0.,0.);
		do {
			//Generate x,y
			double x,y,b;
			gxy(x, y, b, nr);
			//cout << " xcre= " << x << " ycre= " << y << endl;
			double phi,rap1,rap2;
			phi=2.*M_PI*nr.rando();
			rap1=-2.+4.*nr.rando();
			rap2=-2.+4.*nr.rando();

			double l1=0.;
			double l2=0.;
			double L_qhat1=0.;
			double L_qhat2=0.;
			double L2_qhat1=0.;
			double L2_qhat2=0.;
			double L3_qhat1=0.;
			double L3_qhat2=0.;
			double temp1=0.;
			double temp2=0.;
			double twotemp1=0.;
			double twotemp2=0.;
			double ehat1=0.;
			double ehat2=0.;

			double r1[4]={x,y,0.,0.};
			double r2[4]={x,y,0.,0.};

			double p1[4]={cos(phi),sin(phi),sinh(rap1),cosh(rap1)};
			double p2[4]={cos(phi+M_PI),sin(phi+M_PI),sinh(rap2),cosh(rap2)};

			double ang1=acos(p1[0]);
			double ang2=acos(p2[0]);

			double tlab=0.;
		
			bool cold1=0;
			bool cold2=0;

			do {
			  
			  for (unsigned a=0; a<4;a++) {
			    r1[a]+=p1[a]/p1[3]*step;
			    r2[a]+=p2[a]/p2[3]*step;
			  }
		  
   			  double tau1=sqrt(r1[3]*r1[3]-r1[2]*r1[2]);
   			  double tau2=sqrt(r2[3]*r2[3]-r2[2]*r2[2]);
			  double stau1=tau1;
			  double stau2=tau2;
			  if (no_expansion) tau1=0.6, tau2=0.6;

			  if (tau1>=0.6) {

		            vector<double> v;	  
			    double vx=gVx(tau1,r1[0],r1[1],rap1);
			    double vy=gVy(tau1,r1[0],r1[1],rap1);
			    double vz=p1[2]/p1[3];
			    double frap=atanh(vz);
			    vx/=cosh(frap);
			    vy/=cosh(frap);
			    v.push_back(vx), v.push_back(vy), v.push_back(vz), v.push_back(1.);

			    double v2=pow(v[0],2.)+pow(v[1],2.)+pow(v[2],2.);
			    double w2=1.;
	   		    double vscalw=v[0]*p1[0]/p1[3]+v[1]*p1[1]/p1[3]+v[2]*p1[2]/p1[3];
			    double lore=1./sqrt(1.-v2);
			    double f_step = step*sqrt(w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw));
			    if (w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw)<0.) { cout << " l1 problem !! " << endl; exit(1); }
			  

			    double temp=gT(tau1,r1[0],r1[1],rap1);
			    if (temp>Tc) {
			    
			      l1+=f_step;
			      
			      double pscalv = p1[3]*v[3]-p1[0]*v[0]-p1[1]*v[1]-p1[2]*v[2];
			      L_qhat1+= f_step * pscalv/p1[3] * pow(temp,3.);
			      L2_qhat1+= 2. * f_step * l1 * pscalv/p1[3] * pow(temp,3.);
			      L3_qhat1+= 3. * f_step * l1 * pscalv/p1[3] * pow(temp,3.);
			      
			      temp1+= f_step * temp;
			      twotemp1+= f_step * temp * temp;

			      ehat1+= f_step * pscalv/p1[3] * pow(temp,2.);
			    }
			    else if (stau1>18.5) cold1=1;

			  }

			  if (tau2>=0.6) {
	
			    vector<double> v;	  
			    double vx=gVx(tau2,r2[0],r2[1],rap2);
			    double vy=gVy(tau2,r2[0],r2[1],rap2);
			    double vz=p2[2]/p2[3];
			    double frap=atanh(vz);
			    vx/=cosh(frap);
			    vy/=cosh(frap);
			    v.push_back(vx), v.push_back(vy), v.push_back(vz), v.push_back(1.);

			    double v2=pow(v[0],2.)+pow(v[1],2.)+pow(v[2],2.);
			    double w2=1.;
	   		    double vscalw=v[0]*p2[0]/p2[3]+v[1]*p2[1]/p2[3]+v[2]*p2[2]/p2[3];
			    double lore=1./sqrt(1.-v2);
			    double f_step = step*sqrt(w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw));
			    if (w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw)<0.) { cout << " l2 problem !! " << endl; exit(1); }

			    double temp=gT(tau2,r2[0],r2[1],rap2);
			    if (temp>Tc) {
			    
			      l2+=f_step;

			      double pscalv = p2[3]*v[3]-p2[0]*v[0]-p2[1]*v[1]-p2[2]*v[2];
			      L_qhat2+= f_step * std::pow(temp,3.) * pscalv/p2[3];
			      L2_qhat2+= 2. * f_step * l2 * std::pow(temp,3.) * pscalv/p2[3];
			      L3_qhat2+= 3. * f_step * l2 * l2 * std::pow(temp,3.) * pscalv/p2[3];

			      temp2+= f_step * temp;
			      twotemp2+= f_step * temp * temp;
			      
			      ehat2+= f_step * pscalv/p2[3] * pow(temp,2.);
			    
			    }
			    else if (stau2>18.5) cold2=1;

			  }

			} while (!cold1 || !cold2);

			//cout << " l1= " << l1 << " l2= " << l2 << endl;

			// Fill histo
			//
		
			qhatfile << L_qhat1/l1 << " " << L2_qhat1*pow(hbarc,-2.) << " " << L3_qhat1*pow(hbarc,-3.) << " " << l1*pow(hbarc,-1.) << " " << temp1/l1 << " " << twotemp1/l1 << " " << ehat1/l1 << endl;
			qhatfile << L_qhat2/l2 << " " << L2_qhat2*pow(hbarc,-2.) << " " << L3_qhat2*pow(hbarc,-3.) << " " << l2*pow(hbarc,-1.) << " " << temp2/l2 << " " << twotemp2/l2 << " " << ehat2/l2 << endl;

			avevals[0]+=L_qhat1/l1+L_qhat2/l2;
			avevals[1]+=L2_qhat1*pow(hbarc,-2.)+L2_qhat2*pow(hbarc,-2.);
			avevals[2]+=L3_qhat1*pow(hbarc,-3.)+L3_qhat2*pow(hbarc,-3.);
			avevals[3]+=l1*pow(hbarc,-1.)+l2*pow(hbarc,-1.);
			avevals[4]+=temp1/l1+temp2/l2;
			avevals[5]+=twotemp1/l1+twotemp2/l2;
			avevals[6]+=ehat1/l1+ehat2/l2;

			count+=1;
			if (count % 1000 == 0) cout << " Event = " << count << endl;
		} while (count<Nev);
	
		totcount+=1;
	} while (totcount<N);

	qhatfile.close();
	
        for (unsigned a=0; a<7; a++) {
	  aveqhatfile << avevals[a]/2/Nev << " ";
        }
	aveqhatfile << endl;
        aveqhatfile.close();

	return 0;
}
