#pragma once
#include <vector>
#include <string>
#include "ASAM.h"


class CauchyASAM :	public ASAM
{
	double denchy;
	std::vector<double> xmdchy;  /* miscel. parameters */
	std::vector<double> ymdchy;  /* miscel. parameters */
	void Init();

public:
	CauchyASAM(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, int, int, int);
	~CauchyASAM(void){};

	/**** Cauchy ASAM Functions ****/
	void Learn();
	double SAM(double, double);		/* Cauchy SAM */
	double PROB_J(double, double, int);
};

CauchyASAM::CauchyASAM(const std::vector<double>& xvals, const std::vector<double>& yvals, const std::vector<double>& fxyvals, int _numpat, int _numsam, int _numdes):ASAM(xvals, yvals, fxyvals, _numpat, _numsam, _numdes){
	type="Cauchy";
	Init();
	summarizeData();
}



void CauchyASAM::Init(){
	ASAM::Init();
	disps.resize(NUMPAT,2); disps.fill(2*base );      /* "width" of if-part set function */

	/******** Cauchy SAM's parameters ********/
	xmdchy.assign(NUMPAT, 0.0); /* miscel. parameters */
	ymdchy.assign(NUMPAT, 0.0); /* miscel. parameters */	

	/******** Cauchy SAM Initial Learning Rates ********/
	mu_d   = 1.5;
	mu_m   = 1.5;
	mu_cen = 2;
	mu_V   = 1e-7;
}

double CauchyASAM::SAM(double xin, double yin){
	int i;	double av,num, cx, cy;
	denchy = 0.0;
	num = 0.0;

	for (i=0;i<NUMPAT;i++){
		xmdchy[i] = (xin - means(i,0))/disps(i,0);
		ymdchy[i] = (yin - means(i,1))/disps(i,1);
		cx = 1.0/(1.0 + xmdchy[i]*xmdchy[i]);
		cy = 1.0/(1.0 + ymdchy[i]*ymdchy[i]);

		a[i] = cx*cy;
		av = a[i]*volumes[i];
		denchy = denchy+av;
		num = num+av*centroids[i];
	}
	if (denchy != 0.0)  return num/denchy;
	else                return 0.0;
}

double CauchyASAM::PROB_J(double xin, double yin, int j){
	int i;	double av,num, cx, cy;
	denchy = 0.0;
	num = 0.0;

	for (i=0;i<NUMPAT;i++){
		xmdchy[i] = (xin - means(i,0))/disps(i,0);
		ymdchy[i] = (yin - means(i,1))/disps(i,1);
		cx = 1.0/(1.0 + xmdchy[i]*xmdchy[i]);
		cy = 1.0/(1.0 + ymdchy[i]*ymdchy[i]);

		a[i] = cx*cy;
		av = a[i]*volumes[i];
		denchy = denchy+av;
		if (i == j){
			num = av;
		}
	}
	if (denchy != 0.0)  return num/denchy;
	else                return 0.0;
}


void CauchyASAM::Learn(){
	int i,idx;
	double fuzoutchy;
	double cenerr, a_den, pix;
	double dEdF, dEdFa_den;
	double dEdmx, dEddx, dEdmy, dEddy, cenerr_dx, cenerr_dy;
	double xmd, ymd, cx, cy;	

	for (idx=0;idx<NUMSAM;idx++) {
		//Work on training set
		fuzoutchy = SAM(x[idx], y[idx]);
		dEdF = -(sample[idx]-fuzoutchy);
		if (denchy != 0.0) {
			for (i=0;i<NUMPAT;i++) {
				a_den = a[i]/denchy;
				pix = a_den*volumes[i];
				cenerr = centroids[i]-fuzoutchy;

				xmd = (x[idx]-means(i,0))/disps(i,0); ymd = (y[idx]-means(i,1))/disps(i,1);
				cx = 1.0/(1.0 + xmd*xmd);
				cy = 1.0/(1.0 + ymd*ymd);

				cenerr_dx = cenerr/disps(i,0);
				cenerr_dy = cenerr/disps(i,1);
				dEdmx = 2*dEdF*pix*cenerr_dx*xmd*cx*a[i];
				dEdmy = 2*dEdF*pix*cenerr_dy*ymd*cy*a[i];
				dEddx = dEdmx*xmd;
				dEddy = dEdmy*ymd;

				dEdFa_den = dEdF*a_den;

				volumes[i] -= mu_V*dEdFa_den*cenerr/adaptCounter; //volumes[i] = abs(volumes[i]);
				centroids[i] -= mu_cen*dEdF*volumes[i]*a_den/adaptCounter;

				disps(i,0) -= mu_d*dEddx/adaptCounter; disps(i,1) -= mu_d*dEddy/adaptCounter;
				means(i,0) -= mu_m*dEdmx/adaptCounter; means(i,1) -= mu_m*dEdmy/adaptCounter;

				disps(i,0) = abs(disps(i,0)); disps(i,1) = abs(disps(i,1));
			}
		}
	}
}
