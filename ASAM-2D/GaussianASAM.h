#pragma once
#include <vector>
#include <string>
#include "ASAM.h"


class GaussianASAM : public ASAM{

	double dengs;
	std::vector<double> xmdgs;  /* miscel. parameters */
	std::vector<double> ymdgs;  /* miscel. parameters */
	void Init();

public:
	//GaussianASAM(void);
	GaussianASAM(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, int, int, int);
	~GaussianASAM(void){}

	/**** Gaussian ASAM Functions ****/
	void Learn();
	double SAM(double xin, double yin);		/* Gaussian SAM */
	double PROB_J(double xin, double yin, int j);
};


GaussianASAM::GaussianASAM(const std::vector<double>& xvals, const std::vector<double>& yvals, const std::vector<double>& fxyvals, int _numpat, int _numsam, int _numdes):ASAM(xvals, yvals, fxyvals, _numpat, _numsam, _numdes){
	type="Gauss";
	Init();
	summarizeData();
}

//GaussianASAM::~GaussianASAM(void){}

void GaussianASAM::Init(){
	ASAM::Init();
	disps.resize(NUMPAT,2); disps.fill(base);      /* "width" of if-part set function */

	/******** Gaussian SAM's parameters ********/
	xmdgs.assign(NUMPAT, 0.0);  /* miscel. parameters */
	ymdgs.assign(NUMPAT, 0.0);  /* miscel. parameters */	

	/******** Gaussian SAM Initial Learning Rates ********/
	mu_d   = 0.9;
	mu_m   = 0.9;
	mu_cen = 0.1;
	mu_V   = 1e-7;
}

double GaussianASAM::SAM(double xin, double yin){
	int i; double av;
	double num=0.0;
	dengs = 0.0;

	for (i=0;i<NUMPAT;i++){			  //foreach fuzzy rule...
		xmdgs[i] = (xin-means(i,0))/disps(i,0);    //calc intermediate centered gaussian xvar
		ymdgs[i] = (yin-means(i,1))/disps(i,1);    //calc intermediate centered gaussian yvar
		a[i] = exp(-(xmdgs[i]*xmdgs[i]) - (ymdgs[i]*ymdgs[i])); //calc Gaussian fit value
		av = a[i]*volumes[i];               //calc fit-scaled volume of then-part
		dengs = dengs+av;				  //calc denominator of SAM
		num = num+av*centroids[i];			  //calc numerator of SAM
	}
	if (dengs != 0.0)  return num/dengs;
	else               return 0.0;
}

double GaussianASAM::PROB_J(double xin, double yin, int j){
	int i; double av;
	double num=0.0;
	dengs = 0.0;

	for (i=0;i<NUMPAT;i++){			  //foreach fuzzy rule...
		xmdgs[i] = (xin-means(i,0))/disps(i,0);    //calc intermediate centered gaussian xvar
		ymdgs[i] = (yin-means(i,1))/disps(i,1);    //calc intermediate centered gaussian yvar
		a[i] = exp(-(xmdgs[i]*xmdgs[i]) - (ymdgs[i]*ymdgs[i])); //calc Gaussian fit value
		av = a[i]*volumes[i];               //calc fit-scaled volume of then-part
		dengs = dengs+av;				  //calc denominator of SAM
		if ( i == j){
		num = av;			  //calc numerator of SAM
		}
	}
	if (dengs != 0.0)  return num/dengs;
	else               return 0.0;

}


void GaussianASAM::Learn(){ // New variable notation
	int i=0, idx=0;
	double fuzoutgs /*mu_m, mu_d, mu_cen, mu_V,*/;
	double cenerr, a_den, pix;
	double dEdF, dEdFa_den;
	double dEdmx, dEdmy;
	double dEddx, dEddy;

	for (idx=0; idx<NUMSAM; idx++) {
		//Work on training set
		fuzoutgs = SAM(x[idx], y[idx]);	//calc current F[x] 
		dEdF = -(sample[idx] - fuzoutgs);	//calc -(f[x]-F[x]); dEdF = -epsilon
		if (dengs != 0.0) {				//dengs: SAM denom. shared from GaussianSAM
			for (i=0;i<NUMPAT;i++) {
				a_den = a[i]/dengs;
				dEdFa_den = dEdF*a_den;
				pix = a_den*volumes[i];
				cenerr = centroids[i]-fuzoutgs;

				dEdmx = dEdF*pix*cenerr*xmdgs[i]/disps(i,0);
				dEddx = dEdmx*xmdgs[i];

				dEdmy = dEdF*pix*cenerr*ymdgs[i]/disps(i,1);
				dEddy = dEdmy*ymdgs[i];

				volumes[i] -= mu_V*dEdFa_den*cenerr/adaptCounter; //volumes[i] = abs(volumes[i]);
				centroids[i] -= mu_cen*dEdF*a_den*volumes[i]/adaptCounter;   //pix=a_den*volumes[i];

				disps(i,0) -= mu_d*dEddx/adaptCounter; disps(i,1) -= mu_d*dEddy/adaptCounter;
				means(i,0) -= mu_m*dEdmx/adaptCounter; means(i,1) -= mu_m*dEdmy/adaptCounter;

				disps(i,0) = abs(disps(i,0)) ; disps(i,1) = abs(disps(i,1));
			}
		}
	}
	ASAM::Learn();
}
