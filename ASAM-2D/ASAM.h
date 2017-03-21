#pragma once
#include <vector>
#include <string>
#include "Eigen/Core"
#include <iostream>
#include <fstream>

using namespace std;


typedef Eigen::MatrixXd matrix;

/* Non-member utility functions */
std::vector<double> eigvec2stdvec(const Eigen::VectorXd  & vec){
	std::vector<double> fin;
	for(int i = 0; i < vec.size(); i++)	fin.push_back(vec(i));
	return fin;
}

void filefail(std::string out){
	cerr << "Could not open output file: " << out << endl ;
	system("PAUSE");
	exit(1);
}

class ASAM{
	//private:

protected:
	std::string type; //= "(undefined)";
	int NUMPAT;    /* number of fuzzy sets on input or output side*/
	int NUMSAM;    /* use NUMSAM number of points on training grid */
	int NUMDES;    /* use NUMDES number of points on testing grid */
	unsigned int adaptCounter;  /* Adaptation step counter */
	
	double MINX;      /* lower limit for x axis  */
	double MINY;      /* lower limit for y-axis  */
	double MAXX;      /* upper limit for x axis  */
	double MAXY;      /* upper limit for y-axis  */
	double Minf;      /* lower limit for f samples */
	double Maxf;      /* upper limit for f samples */
	std::vector<double> x;      /* x values for training */
	std::vector<double> y;      /* y values for training */
	std::vector<double> xtest;  /* x values for testing */
	std::vector<double> ytest;  /* y values for testing */
	std::vector<double> sample; /* f(x,y) values (little f in paper) for training */
	std::vector<double> des;    /* f(x,y) values (little f in paper) for testing */
	double base;	/* generic nonzero param to use in dispersion initial assignments */

	std::vector<double> centroids;    /* then-part set centroid in Gaussian SAM */
	std::vector<double> volumes;      /* volume (area) of then-part set in Gaussian SAM */
	std::vector<double> a;			  /* set function values */

	matrix means;				/* "location" of gaussian if-part set function */
	matrix disps;				/* "width" of gaussian if-part set function */
	std::vector<double> F;			/* SAM function approximation */
	std::vector<double> Var;		/* SAM function approximation variance*/

	double mu_m, mu_d, mu_cen, mu_V;  /* Adaptive Learning Rates */

	virtual void Init();							/* ASAM Initializer */
	virtual void Learn();							/* ASAM Adaptation */
	virtual double SAM(double xin, double yin)=0;	/* SAM Approximator */
	virtual double PROB_J(double xin, double yin , int j)=0;
	//virtual double VAR(int j);
	//virtual double CENT(int j);

public:
	//ASAM(void); // NOTE!!! No Copy- or Default- Constructor defined
	ASAM(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, int, int, int);
	virtual ~ASAM(void);
	double Approx();								/* SAM MSE calculator */
	void summarizeData();
	void WriteParams(std::string out);
	void WriteEpoch(int epoch);
	void resetCounter();
};

//ASAM::ASAM(){}

ASAM::ASAM(const std::vector<double>& xvals, const std::vector<double>& yvals, const std::vector<double>& fxyvals, int _numpat, int _numsam, int _numdes) : NUMPAT(_numpat),NUMSAM(_numsam),NUMDES(_numdes),type("(undefined)")
{
	if ( (xvals.size() != fxyvals.size()) || (yvals.size() != fxyvals.size()) ){
		std::cerr << "(x,y) <-> f(x,y) Mismatch!!!";
		exit(1);
	}

	adaptCounter = 1;
	xtest.reserve(NUMDES); ytest.reserve(NUMDES); des.reserve(NUMDES);

	float stride = float(NUMDES)/float(NUMSAM) ;  //float stride => non-uniform sampling
	xtest = xvals;  /* x values for testing */
	ytest = yvals;  /* x values for testing */
	des = fxyvals;    
	/* des is f testing matrix.*/

	/* (x,y) values for training */
	/* Selection indices: makes sure float issue don't cause (x,y)->sample mismatches */
	std::vector<unsigned int> index;  /* store linear index of training samples */
	unsigned int i; index.clear();
	for(float ix = 0; ix < xtest.size(); ix+=stride) index.push_back(int(ix));
	x.reserve(index.size()); y.reserve(index.size()); sample.reserve(index.size());

	for(i=0; i < index.size(); i++) x.push_back(xvals[ index[i] ]);
	for(i=0; i < index.size(); i++) y.push_back(yvals[ index[i] ]);

	/* f(x,y) values (little f in paper) for training */
	for(i=0; i < index.size(); i++) sample.push_back(fxyvals[ index[i] ]);
	NUMSAM = index.size();    /* use NUMSAM number of points on training grid */

	MINX = *(std::min_element(xtest.begin(), xtest.end()));      /* lower limit for x-axis  */
	MAXX = *(std::max_element(xtest.begin(), xtest.end()));      /* upper limit for x-axis  */

	MINY = *(std::min_element(ytest.begin(), ytest.end()));      /* lower limit for y-axis  */
	MAXY = *(std::max_element(ytest.begin(), ytest.end()));      /* upper limit for y-axis  */
	base = std::max((double)(MAXX-MINX)/(double)(NUMPAT-1), (double)(MAXY-MINY)/(double)(NUMPAT-1));
	Maxf = *(std::max_element(sample.begin(), sample.end()));
	Minf = *(std::min_element(sample.begin(), sample.end()));
}

ASAM::~ASAM(){
	/*delete x, y, xtest, ytest, index, sample, des, centroids, volumes, a, means, disps, F;*/
}

void ASAM::Init(){
	/******** SAM parameters ********/
	volumes.assign(NUMPAT, 1.0);      /* volume (area) of then-part set in SAM */
	centroids.assign(NUMPAT, 0.0);    /* then-part set centroid in SAM */
	a.assign(NUMPAT, 0.0);		/* set function values */	

	means.resize(NUMPAT,2); //means.fill(0);		/* "location" of if-part set function */
	means.setRandom();	means = (means.array() + 1)*0.5;
	means.col(0) = (means.col(0).array())*(MAXX-MINX)+MINX;
	means.col(1) = (means.col(1).array())*(MAXY-MINY)+MINY;

	disps.resize(NUMPAT,2); disps.fill(base);      /* "width" of if-part set function */
	F.assign(NUMDES, 0.0);   /* SAM approximator values */
	Var.assign(NUMDES, 0.0);   /* SAM approximator values variance */

	Eigen::VectorXd rndcen(NUMPAT);
	rndcen.setLinSpaced(NUMPAT, Minf, Maxf ) ;

	centroids = eigvec2stdvec(rndcen); 
}

void ASAM::Learn(){
	adaptCounter++;
}

double ASAM::Approx(){
	int ind; double err, Fx, var, probj, val,cj, sumerr = 0;
	for(ind=0; ind < NUMDES ; ind++){
		Fx = SAM(xtest[ind], ytest[ind]);
		F[ind] = Fx;		
		err = des[ind] - F[ind];
		sumerr += err*err;
		val = 0.0;
		for (int j = 0; j < NUMPAT; j++){
			probj = PROB_J(xtest[ind], ytest[ind], j);
			var = (double)(Maxf-Minf)/(double)(NUMPAT-1); 
			var = (var/2.0)*(var/2.0);
			cj = centroids[j];
			val = val + (probj*var) + (probj*(cj - Fx)*(cj - Fx));
		}
		Var[ind] = val;
	}
	return sumerr/(double)NUMDES;
}

void ASAM::summarizeData(){
	std::cout << std::endl <<"SAM Type: "<< type << std::endl ;
	std::cout <<  "Testing data set size: "<< des.size() << std::endl ;
	std::cout << "Training data set size: "<< sample.size() << std::endl ;
	std::cout << "Number of Rules: "<< NUMPAT << std::endl ;
	std::cout << "Range of (x,y): (" 
		<< MINX << ", " << MAXX << ") x (" 
		<< MINY << ", " << MAXY << ")" << std::endl;
}

void ASAM::WriteParams(std::string out){
	int i;
	ofstream ofp(out.data(), ios::out);  
	if (ofp.fail()) filefail(out);

	/* Output format: m_x  m_y  d_x  d_y  c  V*/

	for (i=0;i<NUMPAT;i++)
		ofp << means(i,0) << "   " << means(i,1) << "   " << 
		disps(i,0) << "   " << disps(i,1) << "   " << 
		centroids[i]<< "   " << volumes[i] <<std::endl;
	ofp.close();
}

void ASAM::WriteEpoch(int epoch){
	std::ofstream ofp; 
	std::ofstream ofp1;

	std::ostringstream s; //to convert #s into string
	std::ostringstream s1;

	s <<"./" << this->type << "/" << "fuzzyF" << this->type << "-" << epoch << ".dat";
	s1 <<"./" << this->type << "/" << "fuzzyV" << this->type << "-" << epoch << ".dat";

	WriteParams("./" + this->type + "/" + "Parameters.par");

	ofp.open(s.str().c_str(), ios::out); if (ofp.fail()) filefail(s.str());
	ofp1.open(s1.str().c_str(), ios::out); if (ofp1.fail()) filefail(s1.str());

	ofp << "  x   " << "\t \t" << "   y   " << "\t \t" << "  f(x,y)   " << endl;
	ofp1 << "  x   " << "\t \t" << "   y   " << "\t \t" << "V[Z|X= x, Y= y]" << endl;
	for (int k=0;k<NUMDES;k++){ 
		ofp << xtest[k] << "\t" << ytest[k] << "\t"  << F[k] << std::endl;
		ofp1 << xtest[k] << "\t" << ytest[k] << "\t"  << Var[k] << std::endl;
	}
	ofp.close();
	ofp1.close();
}

void ASAM::resetCounter(){adaptCounter=1;}
