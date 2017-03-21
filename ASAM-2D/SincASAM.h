#pragma once
#include <vector>
#include <string>
#include "ASAM.h"


class SincASAM : public ASAM{

	double densinc;
	std::vector<double> xmdsinc;  /* miscel. parameters */
	std::vector<double> ymdsinc;  /* miscel. parameters */
	void Init();

public:
	//SincASAM(void);
	SincASAM(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, int, int, int);
	~SincASAM(void){}

	/****Sinc ASAM Functions****************/
	void Learn();
	double SAM(double xin, double yin);		/* Sinc SAM */
	double PROB_J(double xin, double yin, int j);

};

SincASAM::SincASAM(const std::vector<double>& xvals, const std::vector<double>& yvals, const std::vector<double>& fxyvals, int _numpat, int _numsam, int _numdes):ASAM(xvals, yvals, fxyvals, _numpat, _numsam, _numdes){
	type="Sinc";
	Init();
	summarizeData();
}

//SincASAM::~SincASAM(void){}

void SincASAM::Init(){
	ASAM::Init();
	disps.resize(NUMPAT,2); disps.fill(2*base/3.1416 );      /* "width" of if-part set function */

	/******** Sinc SAM's parameters ********/
	xmdsinc.assign(NUMPAT, 0.0);  /* miscel. parameters */
	ymdsinc.assign(NUMPAT, 0.0);  /* miscel. parameters */
	densinc=0;

	/******** Sinc SAM Initial Learning Rates ********/
	mu_cen = 5e-2;
	mu_V   = 1e-7;
	mu_d   = 0.91;
	mu_m   = 0.91;
}

double SincASAM::SAM(double xin, double yin){
	int i;
	double av=0, num=0, sx, sy;
	densinc = 0;

	for (i=0;i<NUMPAT;i++){
		xmdsinc[i] = (xin - means(i,0))/disps(i,0);
		ymdsinc[i] = (yin - means(i,1))/disps(i,1);

		if (xmdsinc[i] == 0.0) sx = 1.0;
		else sx = sin(xmdsinc[i])/xmdsinc[i]; 
		if(ymdsinc[i] == 0.0) sy = 1.0;
		else sy = sin(ymdsinc[i])/ymdsinc[i];

		a[i] = sx*sy;
		av = a[i]*volumes[i];
		densinc = densinc+av;
		num = num+av*centroids[i];
	}
	if (densinc != 0.0)
		return num/densinc;
	else
		return 0.0;
}


double SincASAM::PROB_J(double xin, double yin, int j){
	int i;
	double av=0, num=0, sx, sy;
	densinc = 0;

	for (i=0;i<NUMPAT;i++){
		xmdsinc[i] = (xin - means(i,0))/disps(i,0);
		ymdsinc[i] = (yin - means(i,1))/disps(i,1);

		if (xmdsinc[i] == 0.0) sx = 1.0;
		else sx = sin(xmdsinc[i])/xmdsinc[i]; 
		if(ymdsinc[i] == 0.0) sy = 1.0;
		else sy = sin(ymdsinc[i])/ymdsinc[i];

		a[i] = sx*sy;
		av = a[i]*volumes[i];
		densinc = densinc+av;

		if (i == j){
		num = av;
		}

	}
	if (densinc != 0.0)
		return num/densinc;
	else
		return 0.0;
}


void SincASAM::Learn(){
	int i, idx;
	double fuzoutsinc;
	double cenerr, cenerr_den, cenerrV_den, cosxmd, cosymd, acosxmdcenerrV_den, acosymdcenerrV_den;
	double dEddx, dEddy, dEdmx, dEdmy, dEdF, dEdFacosblablax, dEdFacosblablay;
	double sx, sy;

	for (idx=0; idx<NUMSAM; idx++) {
		//Work on training set
		fuzoutsinc = SAM(x[idx], y[idx]);	//calc current F[x] 
		dEdF = -(sample[idx] - fuzoutsinc);
		for ( i=0; i<NUMPAT; i++)  {
			if (a[i] != 0.0)  {
				cenerr = centroids[i]-fuzoutsinc;
				cenerr_den = cenerr/densinc;
				cenerrV_den = cenerr_den*volumes[i];

				if (xmdsinc[i] == 0.0) sx = 1.0; else sx = sin(xmdsinc[i])/xmdsinc[i]; 
				if(ymdsinc[i] == 0.0) sy = 1.0;	else sy = sin(ymdsinc[i])/ymdsinc[i];
				// sx = sin(xmdsinc[i])/xmdsinc[i]; sy = sin(ymdsinc[i])/ymdsinc[i];
				cosxmd = cos(xmdsinc[i]); cosymd = cos(ymdsinc[i]); 
				acosxmdcenerrV_den = (a[i]-sy*cosxmd)*cenerrV_den; //Note add'l sy
				acosymdcenerrV_den = (a[i]-sx*cosymd)*cenerrV_den; //Note add'l sx
				dEdFacosblablax = dEdF*acosxmdcenerrV_den;
				dEdFacosblablay = dEdF*acosymdcenerrV_den;
				if (xmdsinc[i] != 0.0)
					dEdmx = dEdFacosblablax/(x[idx]-means(i,0));
				else dEdmx = 0.0;
				if (ymdsinc[i] != 0.0)
					dEdmy = dEdFacosblablay/(y[idx]-means(i,1));
				else dEdmy = 0.0;
				dEddx = dEdFacosblablax/disps(i,0);
				dEddy = dEdFacosblablay/disps(i,1);


				/* Modified Laws */
				centroids[i] -= mu_cen*dEdF*volumes[i]*a[i]/(densinc*adaptCounter); 
				volumes[i] -= mu_V*dEdF*a[i]*cenerr_den/adaptCounter; //volumes[i] = abs(volumes[i]);

				means(i,0) -= mu_m*dEdmx/adaptCounter; disps(i,0) -= mu_d*dEddx/adaptCounter;
				means(i,1) -= mu_m*dEdmy/adaptCounter; disps(i,1) -= mu_d*dEddy/adaptCounter;

				disps(i,0) = abs(disps(i,0)); disps(i,1) = abs(disps(i,1));
			}
		}
	}
	ASAM::Learn();
}
