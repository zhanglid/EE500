#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

#include "SAMUtils.h"

/********Variables, Arrays Defaults*********************/
int NUMPAT = 10;     /* number of fuzzy sets on input or output side*/
int NUMSAM = 200;    /* use NUMSAM number of data pairs for training */
int NUMDES = 400;    /* use NUMDES number of data pairs for testing */
static unsigned int adaptCounter = 1;	/*Adaptation Step Counter*/
double MINX;      /* lower limit for x axis  */
double MAXX;      /* upper limit for x axis  */

double MINY;      /* lower limit for y axis  */
double MAXY;      /* upper limit for y axis  */
vector<double> x;      /* x values for training */
vector<double> sample; /* f(x) values (little f in paper) for training */
vector<double> xtest;  /* x values for testing */
vector<double> des;    /* f(x) values (little f in paper) for testing */
string name[6] = {"Gauss", "Cauchy", "Tanh", "Laplace", "Tri", "Sinc"};
//enum Fitfxn {gauss, cauchy, tanhyp, laplace, triangle, sinc};

/******** Gaussian SAM's parameters ********/
vector<double> mgs;      /* "location" of gaussian if-part set function */
vector<double> dgs;      /* "width" of gaussian if-part set function */
vector<double> cengs;    /* then-part set centroid in Gaussian SAM */
vector<double> Vgs;      /* volume (area) of then-part set in Gaussian SAM */
vector<double> ags;      /* gaussian set values */
double dengs;
vector<double> xmdgs;  /* miscel. parameters */
vector<double> Fgss;    /* Gaussian SAM function approximation */
vector<double> Var_gss;    /* Gaussian SAM function approximation variance */
/*******************************************/
/********  Cauchy SAM's parameters  ********/
vector<double> mchy;     /* "location" of cauchy if-part set function */
vector<double> dchy;     /* "width" of cauchy if-part set function */
vector<double> cenchy;   /* then-part set centroid in Cauchy SAM*/
vector<double> Vchy;     /* volume (area) of then-part set in Cauchy SAM */
vector<double> achy;     /* cauchy set values */
double denchy;
vector<double> xmdchy; /* miscel. parameters */
vector<double> Fchy;    /* Cauchy SAM function approximation */
vector<double> Var_chy;    /* Cauchy SAM function approximation variance*/
/******************************************/
/********   Tanh SAM's parameters  ********/
vector<double> mtanh;    /* "location" of tanh if-part set function */
vector<double> dtanh;    /* parameters in tanh set function */
vector<double> centanh;  /* then-part set centroid in tanh SAM*/
vector<double> Vtanh;    /* volume (area) of then-part set in tanh SAM*/
vector<double> a_tanh;   /* tanh set values */
vector<double> xmdtanh, tanhxmd;      /* miscel. parameters */
vector<double> Ftanh;    /* Hyperbolic tangent SAM function approximation */
vector<double> Var_tanh;    /* Hyperbolic tangent SAM function approximation variance*/
double dentanh;			 /*Tanh SAM denominator*/

/******************************************/
/******** Laplace SAM's parameters ********/
vector<double> mlp;     /* "location" of laplace if-part set function */
vector<double> dlp;     /* "width" of laplace if-part set function */
vector<double> cenlp;   /* then-part set centroid in laplace SAM */
vector<double> Vlp;     /* volume (area) of then-part set in Laplace SAM */
vector<double> alp;    /* laplace set values */
double denlp;
vector<double> fxmdlp;  /* miscel. parameters */
vector<double> Flapl;    /* Laplace SAM function approximation */
vector<double> Var_lapl;    /* Laplace SAM function approximation variance*/
/******************************************/
/******** Triangle SAM's parameters ********/
vector<double> mtri;    /* "location" of triangle set function */
vector<double> dtri;    /* "width" of triangle set function */
vector<double> centri;  /* then-part set centroid in Triangle SAM */
vector<double> Vtri;    /* volume (area) of then-parrt set in Triangle SAM */
vector<double> atri;    /* set function */
double dentri; /* miscel. parameters */
vector<double> Ftri;    /* Triangle SAM function approximation */
vector<double> Var_tri;    /* Triangle SAM function approximation variance*/
/******************************************/
/********** Sinc SAM's parameters **********/
vector<double> msinc;    /* "location" of sinc if-part set function */
vector<double> dsinc;    /* "width" of sinc if-part set function */
vector<double> censinc;  /* then-part set centroid in sinc SAM */
vector<double> Vsinc;    /* volume (area) of then-part set in sinc SAM */
vector<double> asinc;    /* sinc set values */
double densinc;
vector<double> xmdsinc; /* miscel. parameters */
vector<double> Fsinc;    /* Sinc SAM function approximation */
vector<double> Var_sinc;    /* Sinc SAM function approximation variance*/

void InitializeAll(int _numpat, int _numsam, int _numdes)
{
	::NUMPAT = _numpat;     /* number of fuzzy sets on input or output side*/
	::NUMSAM = _numsam;    /* use NUMSAM number of data pairs for training */
	::NUMDES = _numdes;    /* use NUMDES number of data pairs for testing */

	/******** Gaussian SAM's parameters ********/
	::mgs.assign(NUMPAT, 0.0);      /* "location" of gaussian if-part set function */
	::dgs.assign(NUMPAT, 0.0);      /* "width" of gaussian if-part set function */
	::cengs.assign(NUMPAT, 0.0);    /* then-part set centroid in Gaussian SAM */
	::Vgs.assign(NUMPAT, 0.0);      /* volume (area) of then-part set in Gaussian SAM */
	::ags.assign(NUMPAT, 0.0);      /* gaussian set values */
	::xmdgs.assign(NUMPAT, 0.0);  /* miscel. parameters */
	::Fgss.assign(NUMDES, 0.0);    /* Gaussian SAM function approximation */
	::Var_gss.assign(NUMDES, 0.0);  /* Gaussian SAM function approximation variance */
	/*******************************************/
	/********  Cauchy SAM's parameters  ********/
	::mchy.assign(NUMPAT, 0.0);     /* "location" of cauchy if-part set function */
	::dchy.assign(NUMPAT, 0.0);     /* "width" of cauchy if-part set function */
	::cenchy.assign(NUMPAT, 0.0);   /* then-part set centroid in Cauchy SAM*/
	::Vchy.assign(NUMPAT, 0.0);     /* volume (area) of then-part set in Cauchy SAM */
	::achy.assign(NUMPAT, 0.0);     /* cauchy set values */
	::xmdchy.assign(NUMPAT, 0.0); /* miscel. parameters */
	::Fchy.assign(NUMDES, 0.0);    /* Cauchy SAM function approximation */
	::Var_chy.assign(NUMDES, 0.0);  /* Cauchy SAM function approximation variance */
	/******************************************/
	/********   Tanh SAM's parameters  ********/
	::mtanh.assign(NUMPAT, 0.0);    /* "location" of tanh if-part set function */
	::dtanh.assign(NUMPAT, 0.0);    /* parameters in tanh set function */
	::centanh.assign(NUMPAT, 0.0);  /* then-part set centroid in tanh SAM*/
	::Vtanh.assign(NUMPAT, 0.0);    /* volume (area) of then-part set in tanh SAM*/
	::a_tanh.assign(NUMPAT, 0.0);   /* tanh set values */
	::xmdtanh.assign(NUMPAT, 0.0),      /* miscel. parameters */
		::tanhxmd.assign(NUMPAT, 0.0); 
	::Ftanh.assign(NUMDES, 0.0);    /* Hyperbolic tangent SAM function approximation */
	::Var_tanh.assign(NUMDES, 0.0);  /* Tanh SAM function approximation variance */
	/******************************************/
	/******** Laplace SAM's parameters ********/
	::mlp.assign(NUMPAT, 0.0);     /* "location" of laplace if-part set function */
	::dlp.assign(NUMPAT, 0.0);     /* "width" of laplace if-part set function */
	::cenlp.assign(NUMPAT, 0.0);   /* then-part set centroid in laplace SAM */
	::Vlp.assign(NUMPAT, 0.0);     /* volume (area) of then-part set in Laplace SAM */
	::alp.assign(NUMPAT, 0.0);    /* laplace set values */
	::fxmdlp.assign(NUMPAT, 0.0);  /* miscel. parameters */
	::Flapl.assign(NUMDES, 0.0);    /* Laplace SAM function approximation */
	::Var_lapl.assign(NUMDES, 0.0);  /* Laplace SAM function approximation variance */
	/******************************************/
	/******** Triangle SAM's parameters ********/
	::mtri.assign(NUMPAT, 0.0);    /* "location" of triangle set function */
	::dtri.assign(NUMPAT, 0.0);    /* "width" of triangle set function */
	::centri.assign(NUMPAT, 0.0);  /* then-part set centroid in Triangle SAM */
	::Vtri.assign(NUMPAT, 0.0);    /* volume (area) of then-parrt set in Triangle SAM */
	::atri.assign(NUMPAT, 0.0);    /* set function */
	::Ftri.assign(NUMDES, 0.0);    /* Triangle SAM function approximation */
	::Var_tri.assign(NUMDES, 0.0);  /* Triangle SAM function approximation variance */
	/******************************************/
	/********** Sinc SAM's parameters *********/
	::msinc.assign(NUMPAT, 0.0);    /* "location" of sinc if-part set function */
	::dsinc.assign(NUMPAT, 0.0);    /* "width" of sinc if-part set function */
	::censinc.assign(NUMPAT, 0.0);  /* then-part set centroid in sinc SAM */
	::Vsinc.assign(NUMPAT, 0.0);    /* volume (area) of then-part set in sinc SAM */
	::asinc.assign(NUMPAT, 0.0);    /* sinc set values */
	::xmdsinc.assign(NUMPAT, 0.0); /* miscel. parameters */
	::Fsinc.assign(NUMDES, 0.0);    /* Sinc SAM function approximation */
	::Var_sinc.assign(NUMDES, 0.0);  /*Sinc SAM function approximation variance */
	::densinc=0;
	/*******************************************/
}
void InitializeFxn(const vector<double>& xvals, const vector<double>& fxvals){

	if (xvals.size() != fxvals.size()){
		std::cerr << "x <-> f(x) Mismatch!!!";
		exit(1);
	}	
	/********Variables, Arrays*********************/
	::xtest = xvals;  /* x values for testing */
	::des = fxvals;    /* f(x) values (little f in paper) for testing */

	/* x values for training */
	for(unsigned int tx2 = 0; tx2 < xvals.size(); tx2+=2){
		::x.push_back(xvals.at(tx2));}

	/* f(x) values (little f in paper) for training */
	for(unsigned int tfx2 = 0; tfx2 < fxvals.size(); tfx2+=2){
		::sample.push_back(fxvals.at(tfx2));}
	::MINX = xvals.front();      /* lower limit for x axis  */
	::MAXX = xvals.back();      /* upper limit for x axis  */
	int size1 = fxvals.size();
	
	::MINY = vecmin(fxvals);
	::MAXY = vecmax(fxvals);
	
}

void WriteParams(string out, int fxn)
{//enum Fitfxn {gauss, cauchy, tanhyp, laplace, triangle, sinc};
	ofstream ofp;
	int i;
	ofp.open(out.data(), ios::out);  
	if (ofp.fail()) {
		cerr << "Could not open output file\n" ;
		exit(1);
	}

	switch( fxn ) {
	case gauss:
		ofp << " Location"  << "  " << "Dispersion" << "   " << " Centroids" << "   " << "Volume" << endl;
		for (i=0;i<NUMPAT;i++){ofp << mgs[i] << "   " << dgs[i] << "   " << cengs[i]<< "   " << Vgs[i] <<endl;}
		break; 
	case cauchy:	
		ofp << " Location"  << "  " << "Dispersion" << "   " << " Centroids" << "   " << "Volume" << endl;
		for (i=0;i<NUMPAT;i++){ofp << mchy[i] << "   " << dchy[i] << "   " << cenchy[i] << "   " << Vchy[i] <<endl;}
		break;
	case tanhyp:
		ofp << " Location"  << "  " << "Dispersion" << "   " << " Centroids" << "   " << "Volume" << endl;
		for (i=0;i<NUMPAT;i++){
			ofp << mtanh[i] << "   " <<  dtanh[i] << "   "<< centanh[i] << "   " << Vtanh[i] << endl;}
		break;
	case laplace:
		ofp << " Location"  << "  " << "Dispersion" << "   " << " Centroids" << "   " << "Volume" << endl;
		for (i=0;i<NUMPAT;i++){ofp << mlp[i] << "   " << dlp[i] << "   " << cenlp[i] << "   " << Vlp[i] <<endl;}
		break;
	case triangle:
		ofp << " Location"  << "  " << "Dispersion" << "   " << " Centroids" << "   " << "Volume" << endl;
		for (i=0;i<NUMPAT;i++){ofp << mtri[i] << "   " << dtri[i] << "   " << centri[i] << "   " << Vtri[i] <<endl;}
		break;
	case sinc:
		ofp << " Location"  << "  " << "Dispersion" << "   " << " Centroids" << "   " << "Volume" << endl;
		for (i=0;i<NUMPAT;i++){ofp << msinc[i] << "   " << dsinc[i] << "   " << censinc[i] << "   " << Vsinc[i] <<endl;}
		break;
	}
	ofp.close();
}


void WriteEpoch(string basename, int epoch)
{//enum Fitfxn {gauss, cauchy, tanhyp, laplace, triangle, sinc};
	ofstream ofp[6]; 
	ofstream ofp1[6];

	string out[6];
	string vout[6];

	int k;
	int j = 0;
	string bnam, bname;

	for(j=gauss; j<=sinc; j++){ //Modify when Sinc SAM fixed
		ostringstream s;
		s << name[j] << "-" << epoch;
		bnam = basename + s.str() + ".dat";
		bname = "fuzzyV" +s.str() + ".dat";		
		
		out[j] = "./" + name[j] + "/" + bnam ;
		vout[j] = "./" + name[j] + "/" + bname ;
		
		WriteParams("./" + name[j] + "/" + "Parameters.par", j);
		ofp[j].open(out[j].data(), ios::out);
		if (ofp[j].fail()) {
			cerr << "Could not open output file: " << out[j] << endl ;
			system("PAUSE");
			exit(1);
		}

		ofp1[j].open(vout[j].data(), ios::out);
		if (ofp1[j].fail()) {
			cerr << "Could not open output file: " << vout[j] << endl ;
			system("PAUSE");
			exit(1);
		}

		switch( j ) {
		case gauss:
			ofp[j] << " x  "  << "\t" << " F(x)  " << endl;
			ofp1[j] << " x "  << "\t" << "V[Y|X = x] " << endl;
			for (k=0;k<NUMDES;k++) {ofp[j] << xtest[k] << "\t" << Fgss[k] << endl;}
			for (k=0;k<NUMDES;k++) {ofp1[j] << xtest[k] << "\t" << Var_gss[k] << endl;}
			break; 
		case cauchy:
			ofp[j] << " x  "  << "\t" << " F(x)  " << endl;
			ofp1[j] << " x "  << "\t" << "V[Y|X = x] " << endl;
			for (k=0;k<NUMDES;k++) {ofp[j] << xtest[k] << "\t" << Fchy[k] << endl;}
			for (k=0;k<NUMDES;k++) {ofp1[j] << xtest[k] << "\t" << Var_chy[k] << endl;}
			break;
		case tanhyp:
			ofp[j] << " x  "  << "\t" << " F(x)  " << endl;
			ofp1[j] << " x "  << "\t" << "V[Y|X = x] " << endl;
			for (k=0;k<NUMDES;k++) {ofp[j] << xtest[k] << "\t" << Ftanh[k] << endl;}
			for (k=0;k<NUMDES;k++) {ofp1[j] << xtest[k] << "\t" << Var_tanh[k] << endl;}
			break;
		case laplace:
			ofp[j] << " x  "  << "\t" << " F(x)  " << endl;
			ofp1[j] << " x "  << "\t" << "V[Y|X = x] " << endl;
			for (k=0;k<NUMDES;k++) {ofp[j] << xtest[k] << "\t" << Flapl[k] << endl;}
			for (k=0;k<NUMDES;k++) {ofp1[j] << xtest[k] << "\t" << Var_lapl[k] << endl;}
			break;
		case triangle:
			ofp[j] << " x  "  << "\t" << " F(x)  " << endl;
			ofp1[j] << " x "  << "\t" << "V[Y|X = x] " << endl;
			for (k=0;k<NUMDES;k++) {ofp[j] << xtest[k] << "\t" << Ftri[k] << endl;}
			for (k=0;k<NUMDES;k++) {ofp1[j] << xtest[k] << "\t" << Var_tri[k] << endl;}
			break;
		case sinc:
			ofp[j] << " x  "  << "\t" << " F(x)  " << endl;
			ofp1[j] << " x "  << "\t" << "V[Y|X = x] " << endl;
			for (k=0;k<NUMDES;k++) {ofp[j] << xtest[k] << "\t" << Fsinc[k] << endl;}
			for (k=0;k<NUMDES;k++) {ofp1[j] << xtest[k] << "\t" << Var_sinc[k] << endl;}
			break;
		}
		ofp[j].close();
		ofp1[j].close();
	}
}

double vecsum(const vector<double>& vec){
	double sum = 0;
	vector<double>::const_iterator it;
	for ( it = vec.begin(); it < vec.end(); it++ )
		sum += *it;
	return sum;
}

double vecmin(const vector<double>& vec){
	vector<double> svec(vec);
	sort( svec.begin(), svec.end() );
	return svec.front();
}



double vecmax(const vector<double>& vec){
	vector<double> svec(vec);
	sort( svec.begin(), svec.end() );
	return svec.back();
}



double SIGN(double xx)
{
	double ss;
	if (xx > 0.0) ss = 1.0;
	else if (xx < 0.0) ss = -1.0;
	else ss = 0.0;
	return ss;
}

/****Gaussian SAM *************************/

void GaussianSAMlearn()
{
	int i,k;
	double mu_m, mu_d, mu_cen, mu_V, fuzoutgs;
	double cenerr, a_den, pix, cenerr_d;
	double dEdm, dEdd, dEdF, dEdFa_den;

	mu_d   = 1e-2/adaptCounter ;	/* learning rates */
	mu_m   = 1e-2/adaptCounter;	/* modified to an approx harmonic series */
	mu_cen = 1e-2/adaptCounter;	/* sum(mu(t)) ~ 1/t = inf */
	mu_V   = 1e-2/adaptCounter;	/* sum(mu(t)^2) ~ 1/t^2 < inf */

	for (k=0;k<NUMSAM;k++) {			//Work on training set
		fuzoutgs = GaussianSAM(x[k]);	//calc current F[x] 
		dEdF = -(sample[k]-fuzoutgs);	//calc -(f[x]-F[x])
		if (dengs != 0.0) {				//dengs: SAM denom. shared from GaussianSAM
			for (i=0;i<NUMPAT;i++) {		//Update params foreach Rule using error info
				a_den = ags[i]/dengs;
				pix = a_den*Vgs[i];
				cenerr = cengs[i]-fuzoutgs;
				cenerr_d = cenerr/dgs[i];

				dEdm = dEdF*pix*cenerr_d*xmdgs[i];
				dEdd = dEdm*xmdgs[i];

				dEdFa_den = dEdF*a_den;		//Not quite dEdc[i]. missing V[i] factor

				dgs[i] -= mu_d*dEdd;
				mgs[i] -= mu_m*dEdm;
				cengs[i] -= mu_cen*Vgs[i]*dEdFa_den;	//-=mu_c*dEdc[i]
				Vgs[i] -= mu_V*dEdFa_den*cenerr;
			}
		}
	}
}


double GaussianSAM(double in)
{
	int i;
	double av, num=0.0;

	dengs = 0.0;
	for (i=0;i<NUMPAT;i++)  {			  //foreach fuzzy rule...
		xmdgs[i] = (in-mgs[i])/dgs[i];    //calc intermediate centered gaussian var
		ags[i] = exp(-xmdgs[i]*xmdgs[i]); //calc Gaussian fit value
		av = ags[i]*Vgs[i];               //calc fit-scaled volume of then-part
		dengs = dengs+av;				  //calc denominator of SAM
		num = num+av*cengs[i];			  //calc numerator of SAM
	}
	if (dengs != 0.0)  return num/dengs;
	else               return 0.0;
}


double GaussianSAMProb(double in, int j)
{
	int i;
	double av, num=0.0;

	dengs = 0.0;
	for (i=0;i<NUMPAT;i++)  {		  //foreach fuzzy rule...
		xmdgs[i] = (in-mgs[i])/dgs[i];    //calc intermediate centered gaussian var
		ags[i] = exp(-xmdgs[i]*xmdgs[i]); //calc Gaussian fit value
		av = ags[i]*Vgs[i];               //calc fit-scaled volume of then-part
		dengs = dengs+av;		  //calc denominator of SAM
		
		if (i == j){
		    num = av;           //calc numerator of SAM
		}
	}
	
	if (dengs != 0.0)  return num/dengs;
	else               return 0.0;
}


double GaussianSAMVar(int j)
{
	double base = (double)(MAXY-MINY)/(double)(NUMPAT-1);
	return (base/1.52)*(base/1.52);
}


double GaussianSAMapprox()         //returns MSE of current approximation
{
	int j, k;
	double err, sumerr= 0.0;
	double probj, var, varj, Fx;

	for (k=0;k<NUMDES;k++) {
		Fx = GaussianSAM(xtest[k]); //Updates approximation using SAM computed
		Fgss[k] = Fx;		
		err = des[k]-Fgss[k];
		sumerr += err*err;

		// Compute the conditional variance
		var = 0.0;
		for (j = 0; j < NUMPAT; j++){
			probj = GaussianSAMProb(xtest[k], j);
			varj  = GaussianSAMVar(j);
			var = ((probj*varj) + (probj*(cengs[j] - Fx)*(cengs[j] - Fx))) + var;
		}
		Var_gss[k] = var;
	}
	return sumerr/(double)NUMDES;
}


void GaussianSAMinit()
{
	int i;
	double j,step,base;

	step = ((double)(NUMSAM-1))/((double)(NUMPAT-1));
	base = (double)(MAXX-MINX)/(double)(NUMPAT-1);
	i = -1;
	for (j=0.0;j<(double)NUMSAM;j+=step) {
		i += 1;
		cengs[i] = sample[(int)(j+.5)];
		dgs[i] = base/1.52;
		mgs[i] = x[(int)(j+.5)];
		Vgs[i] = 1.0;
	}
	
}


/****Sinc SAM*******************************/
void SincSAMlearn()
{
	int i,k;
	double mu_m, mu_d, mu_cen, mu_V, fuzoutsinc;
	double cenerr, cenerr_den, cenerrV_den, cosxmd, acosxmdcenerrV_den;
	double dEdd, dEdm, dEdF, dEdFacosblabla;

	mu_d   = 1e-5/adaptCounter ;	/* learning rates */
	mu_m   = 1e-5/adaptCounter;	/* modified to an approx harmonic series */
	mu_cen = 1e-5/adaptCounter;	/* sum(mu(t)) ~ 1/t = inf */
	mu_V   = 1e-5/adaptCounter;	/* sum(mu(t)^2) ~ 1/t^2 < inf */

	for (k=0;k<NUMSAM;k++)  {
		fuzoutsinc = SincSAM(x[k]);
		dEdF = -(sample[k]-fuzoutsinc);
		for (i=0;i<NUMPAT;i++)  {
			if (asinc[i] != 0.0)  {
				cenerr = censinc[i]-fuzoutsinc;
				cenerr_den = cenerr/densinc;
				cenerrV_den = cenerr_den*Vsinc[i];
				cosxmd = cos(xmdsinc[i]);
				acosxmdcenerrV_den = (asinc[i]-cosxmd)*cenerrV_den;
				dEdFacosblabla = dEdF*acosxmdcenerrV_den;
				if (xmdsinc[i] != 0.0)
					dEdm = dEdFacosblabla/(x[k]-msinc[i]);
				else dEdm = 0.0;
				dEdd = dEdFacosblabla/dsinc[i];

				dsinc[i] -= mu_d*dEdd;
				msinc[i] -= mu_m*dEdm;
				censinc[i] -= mu_cen*dEdF*Vsinc[i]*asinc[i]/densinc; 
				Vsinc[i] -= mu_V*dEdF*asinc[i]*cenerr_den;
			}
		}
	}
}

double SincSAM(double in)
{
	int i;
	double av=0, num=0;

	densinc = 0;
	for (i=0;i<NUMPAT;i++)  {
		xmdsinc[i] = (in - msinc[i])/dsinc[i];
		if (xmdsinc[i] == 0.0)  asinc[i] = 1.0;
		else  asinc[i] = sin(xmdsinc[i])/xmdsinc[i];

		av = asinc[i]*Vsinc[i];
		densinc = densinc+av;
		num = num+av*censinc[i];
	}
	if (densinc > 0.0)  return num/densinc;
	else if (densinc < 0.0)  {
		std::cout << "densinc = " << densinc << " < 0!!\n";
		return 0.0;
	}
	else                 return 0.0;
}


double SincSAMProb(double in, int j)
{
	int i;
	double av=0, num=0;

	densinc = 0;
	for (i=0;i<NUMPAT;i++)  {
		xmdsinc[i] = (in - msinc[i])/dsinc[i];
		if (xmdsinc[i] == 0.0)  asinc[i] = 1.0;
		else  asinc[i] = sin(xmdsinc[i])/xmdsinc[i];
		av = asinc[i]*Vsinc[i];
		densinc = densinc+av;

		if (i == j){
		    num = av;           //calc numerator of SAM
		}

	}

	if (densinc > 0.0)  return num/densinc;
	else if (densinc < 0.0)  {
		std::cout << "densinc = " << densinc << " < 0!!\n";
		return 0.0;
	}
	else                 return 0.0;
}


double SincSAMVar(int j)
{
	double base = (double)(MAXY-MINY)/(double)(NUMPAT-1);
	return (base/3.1416)*(base/3.1416);
}


void SincSAMinit()
{
	int i;
	double j,STEP,base;
	STEP = ((double)(NUMSAM-1))/((double)(NUMPAT-1));
	base = (double)(MAXX-MINX)/(double)(NUMPAT-1);
	i = -1;
	for (j=0.0;j<(double)NUMSAM;j+=STEP)  {
		i += 1;
		censinc[i] = sample[(int)(j)];
		dsinc[i] = base/3.1416;
		msinc[i] = x[(int)(j)];
		Vsinc[i] = 1.0;
	}
}


double SincSAMapprox()
{

	int j, k;
	double err, sumerr=0.0;
	double probj, var, varj, Fx;

	for (k=0;k<NUMDES;k++) {
		Fx = SincSAM(xtest[k]);
		Fsinc[k] = Fx;
		err = des[k]-Fsinc[k];
		sumerr += err*err;

		// Compute the conditional variance
		var = 0.0;
		for (j = 0; j < NUMPAT; j++){
			probj = SincSAMProb(xtest[k], j);
			varj  = SincSAMVar(j);
			var = ((probj*varj) + (probj*(censinc[j] - Fx)*(censinc[j] - Fx))) + var;
		}
		Var_sinc[k] = var;
	}
	return sumerr/(double)NUMDES;
}


/****Cauchy SAM*******************************/
void CauchySAMlearn()
{
	int i,k;
	double mu_m, mu_d, mu_cen, mu_V, fuzoutchy;
	double cenerr, a_den, pix, cenerr_d;
	double dEdm, dEdd, dEdF, dEdFa_den;

	mu_d   = 2e-2/adaptCounter ;	/* learning rates */
	mu_m   = 8e-3/adaptCounter;	/* modified to an approx harmonic series */
	mu_cen = 2e-1/adaptCounter;	/* sum(mu(t)) ~ 1/t = inf */
	mu_V   = 2e-1/adaptCounter;	/* sum(mu(t)^2) ~ 1/t^2 < inf */

	for (k=0;k<NUMSAM;k++) {
		fuzoutchy = CauchySAM(x[k]);
		dEdF = -(sample[k]-fuzoutchy);
		if (denchy != 0.0) {
			for (i=0;i<NUMPAT;i++) {
				a_den = achy[i]/denchy;
				pix = a_den*Vchy[i];
				cenerr = cenchy[i]-fuzoutchy;
				cenerr_d = cenerr/dchy[i];

				dEdm = dEdF*pix*cenerr_d*xmdchy[i]*achy[i];
				dEdd = dEdm*xmdchy[i];

				dEdFa_den = dEdF*a_den;

				dchy[i] -= mu_d*dEdd;
				mchy[i] -= mu_m*dEdm;
				cenchy[i] -= mu_cen*Vchy[i]*dEdFa_den;
				Vchy[i] -= mu_V*dEdFa_den*cenerr;
			}
		}
	}
}

double CauchySAM(double in)
{
	int i;
	double av,num;

	denchy = 0.0;
	num = 0.0;
	for (i=0;i<NUMPAT;i++)  {
		xmdchy[i] = (in - mchy[i])/dchy[i];
		achy[i] = 1.0/(1.0 + xmdchy[i]*xmdchy[i]);

		av = achy[i]*Vchy[i];
		denchy = denchy+av;
		num = num+av*cenchy[i];
	}
	if (denchy != 0.0)  return num/denchy;
	else                return 0.0;
}


double CauchySAMProb(double in, int j)
{
	int i;
	double av,num;

	denchy = 0.0;
	num = 0.0;
	for (i=0;i<NUMPAT;i++)  {
		xmdchy[i] = (in - mchy[i])/dchy[i];
		achy[i] = 1.0/(1.0 + xmdchy[i]*xmdchy[i]);

		av = achy[i]*Vchy[i];
		denchy = denchy + av;

		if (i == j){
		    num = av;           //calc numerator of SAM
		}
	}

	if (denchy != 0.0)  return num/denchy;
	else                return 0.0;
}


double CauchySAMVar(int j)
{
	double base = (double)(MAXY-MINY)/(double)(NUMPAT-1);
	return (base/3.0)*(base/3.0);
}


void CauchySAMinit()
{
	int i;
	double j,STEP,base;
	STEP = ((double)(NUMSAM-1))/((double)(NUMPAT-1));
	base = (double)(MAXX-MINX)/(double)(NUMPAT-1);
	i = -1;
	for (j=0.0;j<(double)NUMSAM;j+=STEP) {
		i += 1;
		cenchy[i] = sample[(int)(j+.5)];
		dchy[i] = base/3.0;
		mchy[i] = x[(int)(j+.5)];
		Vchy[i] = 1.0;
	}
}


double CauchySAMapprox()
{
	int j, k;
	double err, sumerr= 0.0;
	double probj, var, varj, Fx;

	for (k=0;k<NUMDES;k++) {
		Fx = CauchySAM(xtest[k]);
		Fchy[k] = Fx;
		err = des[k]-Fchy[k];
		sumerr += err*err;

		// Compute the conditional variance
		var = 0.0;
		for (j = 0; j < NUMPAT; j++){
			probj = CauchySAMProb(xtest[k], j);
			varj  = CauchySAMVar(j);
			var = ((probj*varj) + (probj*(cenchy[j] - Fx)*(cenchy[j] - Fx))) + var;
		}
		Var_chy[k] = var;
	}
	return sumerr/(double)NUMDES;
}


/****Tanh SAM********************************/
void TanhSAMlearn()
{
	int i,k;
	double mu_m, mu_d, mu_cen, mu_V, fuzouttanh;
	double a_den, cenerr, cenerr_d, pix;
	double dEdF /*, dFdaj ,dEdFa_den,*/;
	//double dajdd, dajdm;

	mu_d   = 2e-3/adaptCounter ;	/* learning rates */
	mu_m   = 2e-3/adaptCounter;	/* modified to an approx harmonic series */
	mu_cen = 1e-3/adaptCounter;	/* sum(mu(t)) ~ 1/t = inf */
	mu_V   = 1e-3/adaptCounter;	/* sum(mu(t)^2) ~ 1/t^2 < inf */

	for (k=0;k<NUMSAM;k++)  {
		fuzouttanh = TanhSAM(x[k]);
		dEdF = -(sample[k]-fuzouttanh);
		if (dentanh != 0.0)  {
			for (i=0;i<NUMPAT;i++)  {
				if (a_tanh[i] > 0.0) {
					cenerr = centanh[i]-fuzouttanh;
					cenerr_d = cenerr/dtanh[i];
					a_den = a_tanh[i]/dentanh;
					pix = a_den*Vtanh[i];
					//dEdFa_den = dEdF*a_den;

					mtanh[i] -= mu_m*dEdF*pix*cenerr_d*(2-a_tanh[i])*xmdtanh[i];
					dtanh[i] -= mu_d*dEdF*pix*cenerr_d*(2-a_tanh[i])*xmdtanh[i]*xmdtanh[i];
					centanh[i] -= mu_cen*dEdF*pix;
					Vtanh[i] -= mu_V*dEdF*cenerr*a_den;
				}
			}
		}
	}
}


double TanhSAM(double in)
{
	int i;
	double xm, av, num;

	dentanh = 0.0;
	num = 0.0;
	for (i=0;i<NUMPAT;i++) {
		xm = in - mtanh[i];
		xmdtanh[i] = xm/dtanh[i];
		tanhxmd[i] = tanh(-(xmdtanh[i]*xmdtanh[i]) );
		a_tanh[i] = 1 + tanhxmd[i];
		av = a_tanh[i]*Vtanh[i];
		dentanh = dentanh+av;
		num = num+av*centanh[i];
	}
	if (dentanh != 0.0)  return num/dentanh;
	else                 return 0.0;
}


double TanhSAMProb(double in, int j)
{
	int i;
	double xm, av, num;

	dentanh = 0.0;
	num = 0.0;
	for (i=0;i<NUMPAT;i++) {
		xm = in - mtanh[i];
		xmdtanh[i] = xm/dtanh[i];
		tanhxmd[i] = tanh(-(xmdtanh[i]*xmdtanh[i]) );
		a_tanh[i] = 1 + tanhxmd[i];
		av = a_tanh[i]*Vtanh[i];
		dentanh = dentanh+av;

		if ( i == j) {
		    num = av;      // Compute the numerator
		}
	}
	if (dentanh != 0.0)  return num/dentanh;
	else                 return 0.0;
}


double TanhSAMVar(int j)
{
	double base = (double)(MAXY-MINY)/(double)(NUMPAT-1);
	return base*base;
}


void TanhSAMinit()
{
	int i;
	double j,STEP,base;

	
	STEP = ((double)(NUMSAM-1))/((double)(NUMPAT-1));
	base = (double)(MAXX-MINX)/(double)(NUMPAT-1);
	i = -1;
	for (j=0.0;j<(double)NUMSAM;j+=STEP)  {
		i += 1;
		dtanh[i] = base;
		mtanh[i] = x[(int)(j+.5)];
		centanh[i] = sample[(int)(j+.5)];
		Vtanh[i] = 1.0;
	}
}


double TanhSAMapprox()
{	
	int j, k;
	double err, sumerr=0.0;
	double probj, var, varj, Fx;

	for (k=0;k<NUMDES;k++) {
		Fx = TanhSAM(xtest[k]);
		Ftanh[k] = Fx;
		err = des[k]-Ftanh[k];
		sumerr += err*err;

		// Compute the conditional variance
		var = 0.0;
		for (j = 0; j < NUMPAT; j++){
			probj = TanhSAMProb(xtest[k], j);
			varj  = TanhSAMVar(j);
			var = ((probj*varj) + (probj*(centanh[j] - Fx)*(centanh[j] - Fx))) + var;
		}
		Var_tanh[k] = var;
	}
	return sumerr/(double)NUMDES;
}


/****Laplace SAM*****************************/
void LaplaceSAMlearn()
{
	int i,k;
	double mu_m, mu_d, mu_cen, mu_V, fuzoutlp;
	double cenerr, err, xm, dEdm, dEdd, dEdFpicenerr_fd;
	double pj;  /* convex coefficient */

	mu_d   = 2e-2/adaptCounter;	/* learning rates */
	mu_m   = 8e-3/adaptCounter;	/* modified to an approx harmonic series */
	mu_cen = 2e-1/adaptCounter;	/* sum(mu(t)) ~ 1/t = inf */
	mu_V   = 2e-1/adaptCounter;	/* sum(mu(t)^2) ~ 1/t^2 < inf */

	for (k=0;k<NUMSAM;k++)  {
		fuzoutlp = LaplaceSAM(x[k]);
		err = sample[k]-fuzoutlp;
		if (denlp != 0.0)  {
			for (i=0;i<NUMPAT;i++)  {
				xm = x[k]-mlp[i];
				cenerr = cenlp[i] - fuzoutlp;
				pj = Vlp[i]*alp[i]/denlp;
				dEdFpicenerr_fd = -err*pj*cenerr/fabs(dlp[i]);
				dEdm = dEdFpicenerr_fd*SIGN(xm);
				dEdd = dEdFpicenerr_fd*SIGN(dlp[i])*fxmdlp[i];

				mlp[i] -=  mu_m*dEdm;  
				dlp[i] -= mu_d*dEdd;
				Vlp[i] += mu_V*err*alp[i]*cenerr/denlp;
				cenlp[i] += mu_cen*err*pj;
			}
		}
	}
}

double LaplaceSAM(double in)
{
	int i;
	double xm, av, num;

	denlp = 0.0;
	num = 0.0;
	for (i=0;i<NUMPAT;i++) {
		xm = in-mlp[i];
		fxmdlp[i] = fabs(xm/dlp[i]);
		alp[i] = exp(-fxmdlp[i]);
		av = alp[i]*Vlp[i];
		denlp = denlp+av;
		num = num+av*cenlp[i];
	}
	if (denlp != 0.0)  return num/denlp;
	else               return 0.0;
}

double LaplaceSAMProb(double in, int j)
{
	int i;
	double xm, av, num;

	denlp = 0.0;
	num = 0.0;
	for (i=0;i<NUMPAT;i++) {
		xm = in-mlp[i];
		fxmdlp[i] = fabs(xm/dlp[i]);
		alp[i] = exp(-fxmdlp[i]);
		av = alp[i]*Vlp[i];
		denlp = denlp+av;

		if ( i == j) {
		    num = av;      // Compute the numerator
		}
	}
	if (denlp != 0.0)  return num/denlp;
	else               return 0.0;
}


double LaplaceSAMVar(int j)
{
	double base = (double)(MAXY-MINY)/(double)(NUMPAT-1);
	return base/3.0;
}

void LaplaceSAMinit()
{
	int i;
	double j,STEP,base;
	STEP = ((double)(NUMSAM-1))/((double)(NUMPAT-1));
	base = (double)(MAXX-MINX)/(double)(NUMPAT-1);
	i = -1;
	for (j=0.0;j<(double)NUMSAM;j+=STEP) {
		i += 1;
		cenlp[i] = sample[(int)(j+.5)];
		dlp[i] = base/3.0;
		mlp[i] = x[(int)(j+.5)];
		Vlp[i] = 1.0;
	}
}

double LaplaceSAMapprox()
{
	int j, k;
	double err, sumerr=0.0;
	double probj, var, varj, Fx;

	for (k=0;k<NUMDES;k++) {
		Fx = LaplaceSAM(xtest[k]);
		Flapl[k] = Fx;
		err = des[k]-Flapl[k];
		sumerr += err*err;

		// Compute the conditional variance
		var = 0.0;
		for (j = 0; j < NUMPAT; j++){
			probj = LaplaceSAMProb(xtest[k], j);
			varj  = LaplaceSAMVar(j);
			var = ((probj*varj) + (probj*(cenlp[j] - Fx)*(cenlp[j] - Fx))) + var;
		}
		Var_lapl[k] = var;
	}
	return sumerr/(double)NUMDES;
}



/*****Triangle SAM****************************/
void TriangleSAMlearn()  /* Use Gaussian learning laws to tune m and d */
{
	int i,k;
	double mu_d, mu_m, mu_cen, mu_V, fuzouttri;
	double cenerr, err, dEdd, dEdm;
	double pj;  /* convex coefficient */

	mu_d   = 2e-3/adaptCounter ;	/* learning rates */
	mu_m   = 2e-3/adaptCounter;	/* modified to an approx harmonic series */
	mu_cen = 1e-3/adaptCounter;	/* sum(mu(t)) ~ 1/t = inf */
	mu_V   = 1e-3/adaptCounter;	/* sum(mu(t)^2) ~ 1/t^2 < inf */

	for (k=0;k<NUMSAM;k++)  {
		fuzouttri = TriangleSAM(x[k]);
		err = sample[k]-fuzouttri;
		if (dentri != 0.0)  {
			for (i=0;i<NUMPAT;i++)  {
				cenerr = centri[i] - fuzouttri;
				pj = Vtri[i]*atri[i]/dentri;
				dEdm = - err*pj*cenerr*(x[k]-mtri[i])/(dtri[i]*dtri[i]);
				dEdd = dEdm*(x[k]-mtri[i])/dtri[i];
				dtri[i] -= mu_d*dEdd;
				mtri[i] -= mu_m*dEdm;
				Vtri[i] += mu_V*err*atri[i]*cenerr/dentri;
				centri[i] += mu_cen*err*pj;
			}
		}
	}
}

double TriangleSAM(double in)
{
	int i;
	double av, num, fabs_delta;

	dentri = 0.0;
	num = 0.0;

	for (i=0;i<NUMPAT;i++)  {
		fabs_delta = fabs(in - mtri[i]);
		if (fabs_delta <= dtri[i])
			atri[i] = 1.0 - fabs_delta/dtri[i];
		else atri[i] = 0.0;
		av = atri[i]*Vtri[i];
		dentri = dentri+av;
		num = num+av*centri[i];
	}
	if (dentri != 0.0)  return num/dentri;
	else                return 0.0;
}


double TriangleSAMProb(double in, int j)
{
	int i;
	double av, num, fabs_delta;

	dentri = 0.0;
	num = 0.0;

	for (i=0;i<NUMPAT;i++)  {
		fabs_delta = fabs(in - mtri[i]);
		if (fabs_delta <= dtri[i])
			atri[i] = 1.0 - fabs_delta/dtri[i];
		else atri[i] = 0.0;
		av = atri[i]*Vtri[i];
		dentri = dentri+av;
		if (i == j){
		    num = av;
		}
	}

	if (dentri != 0.0)  return num/dentri;
	else                return 0.0;	
}


double TriangleSAMVar(int j)
{
	double base = (double)(MAXY-MINY)/(double)(NUMPAT-1);
	return base;	
}


void TriangleSAMinit()
{
	int i;
	double j,STEP,base;
	STEP = (double)(NUMSAM-1)/(double)(NUMPAT-1);
	base = (double)(MAXX-MINX)/(double)(NUMPAT-1);
	i = -1;
	for (j=0.0;j<(double)NUMSAM;j+=STEP) {
		i += 1;
		centri[i] = sample[(int)(j+.5)];
		dtri[i] = base;
		mtri[i] = x[(int)(j+.5)];
		Vtri[i] = 1.0;
	}
}

double TriangleSAMapprox(){

	int j, k;
	double err, sumerr=0.0;
	double probj, var, varj, Fx;

	for (k=0;k<NUMDES;k++) {
		Fx = TriangleSAM(xtest[k]);
		Ftri[k] = Fx;
		err = des[k]-Ftri[k];
		sumerr += err*err;

		// Compute the conditional variance
		var = 0.0;
		for (j = 0; j < NUMPAT; j++){
			probj = TriangleSAMProb(xtest[k], j);
			varj  = TriangleSAMVar(j);
			var = ((probj*varj) + (probj*(centri[j] - Fx)*(centri[j] - Fx))) + var;
		}
		Var_tri[k] = var;
	}
	return sumerr/(double)NUMDES;
}


/******Aggregate ASAM Calls**********************/
void ASAMsInitialize(){
	adaptCounter = 1;
	GaussianSAMinit();
	SincSAMinit();  
	CauchySAMinit(); 
	TanhSAMinit(); 
	LaplaceSAMinit(); 
	TriangleSAMinit(); 
}

void ASAMsLearn(){
	GaussianSAMlearn();
	SincSAMlearn();  
	CauchySAMlearn(); 
	TanhSAMlearn(); 
	LaplaceSAMlearn(); 
	TriangleSAMlearn();
	adaptCounter++;		// Advance adaptation clock
	/*if (adaptCounter % 1000 == 0) std::cout << "\n Adaptation step # "
	<< adaptCounter << endl ;*/
}

vector<double> ASAMsApprox(){// err vector fill order important!
	vector<double> err;
	//enum Fitfxn {gauss, cauchy, tanhyp, laplace, triangle, sinc};
	err.push_back( GaussianSAMapprox() );
	err.push_back( CauchySAMapprox() ); 
	err.push_back( TanhSAMapprox() ); 
	err.push_back( LaplaceSAMapprox()); 
	err.push_back( TriangleSAMapprox() ); 
	err.push_back( SincSAMapprox() );
	return err;
}
