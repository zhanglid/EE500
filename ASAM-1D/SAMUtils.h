#ifndef SAMUTILS_H
#define SAMUTILS_H

/*******Utility Functions*********************/
enum Fitfxn {gauss, cauchy, tanhyp, laplace, triangle, sinc};
void InitializeAll(int numpat, int numsam, int numdes);
void InitializeFxn(const vector<double>& xvals, const vector<double>& fxvals);
void WriteEpoch(string basename, int epoch);
void WriteParams(string out, int fxn);
double vecsum(const vector<double>& vec);
double vecmin(const vector<double>& vec);
double vecmax(const vector<double>& vec);
double SIGN(double xx);

/********Variables, Arrays Defaults*********************/
extern int NUMPAT;    /* number of fuzzy sets on input or output side*/
extern int NUMSAM;    /* use NUMSAM number of data pairs for training */
extern int NUMDES;    /* use NUMDES number of data pairs for testing */
extern double MINX;      /* lower limit for x axis  */
extern double MAXX;      /* upper limit for x axis  */
extern vector<double> x;      /* x values for training */
extern vector<double> sample; /* f(x) values (little f in paper) for training */
extern vector<double> xtest;  /* x values for testing */
extern vector<double> des;    /* f(x) values (little f in paper) for testing */
extern string name[];

/******** Gaussian SAM's parameters ********/
extern vector<double> mgs;      /* "location" of gaussian if-part set function */
extern vector<double> dgs;      /* "width" of gaussian if-part set function */
extern vector<double> cengs;    /* then-part set centroid in Gaussian SAM */
extern vector<double> Vgs;      /* volume (area) of then-part set in Gaussian SAM */
extern vector<double> ags;      /* gaussian set values */
extern double dengs;
extern vector<double> xmdgs;  /* miscel. parameters */
extern vector<double> Fgss;    /* Gaussian SAM function approximation */
/*******************************************/
/********  Cauchy SAM's parameters  ********/
extern vector<double> mchy;     /* "location" of cauchy if-part set function */
extern vector<double> dchy;     /* "width" of cauchy if-part set function */
extern vector<double> cenchy;   /* then-part set centroid in Cauchy SAM*/
extern vector<double> Vchy;     /* volume (area) of then-part set in Cauchy SAM */
extern vector<double> achy;     /* cauchy set values */
extern double denchy;
extern vector<double> xmdchy; /* miscel. parameters */
extern vector<double> Fchy;    /* Cauchy SAM function approximation */
/******************************************/
/********   Tanh SAM's parameters  ********/
extern vector<double> mtanh;    /* "location" of tanh if-part set function */
extern vector<double> dtanh;    /* parameters in tanh set function */
extern vector<double> centanh;  /* then-part set centroid in tanh SAM*/
extern vector<double> Vtanh;    /* volume (area) of then-part set in tanh SAM*/
extern vector<double> a_tanh;   /* tanh set values */
extern vector<double> xmdtanh, tanhxmd;      /* miscel. parameters */
extern vector<double> Ftanh;    /* Hyperbolic tangent SAM function approximation */
extern double dentanh;			 /*Tanh SAM denominator*/
/******************************************/
/******** Laplace SAM's parameters ********/
extern vector<double> mlp;     /* "location" of laplace if-part set function */
extern vector<double> dlp;     /* "width" of laplace if-part set function */
extern vector<double> cenlp;   /* then-part set centroid in laplace SAM */
extern vector<double> Vlp;     /* volume (area) of then-part set in Laplace SAM */
extern vector<double> alp;    /* laplace set values */
extern double denlp;
extern vector<double> fxmdlp;  /* miscel. parameters */
extern vector<double> Flapl;    /* Laplace SAM function approximation */
/******************************************/
/******** Triangle SAM's parameters ********/
extern vector<double> mtri;    /* "location" of triangle set function */
extern vector<double> dtri;    /* "width" of triangle set function */
extern vector<double> centri;  /* then-part set centroid in Triangle SAM */
extern vector<double> Vtri;    /* volume (area) of then-parrt set in Triangle SAM */
extern vector<double> atri;    /* set function */
extern double dentri; /* miscel. parameters */
extern vector<double> Ftri;    /* Triangle SAM function approximation */
/******************************************/
/********** Sinc SAM's parameters **********/
extern vector<double> msinc;    /* "location" of sinc if-part set function */
extern vector<double> dsinc;    /* "width" of sinc if-part set function */
extern vector<double> censinc;  /* then-part set centroid in sinc SAM */
extern vector<double> Vsinc;    /* volume (area) of then-part set in sinc SAM */
extern vector<double> asinc;    /* sinc set values */
extern double densinc;
extern vector<double> xmdsinc; /* miscel. parameters */
extern vector<double> Fsinc;    /* Sinc SAM function approximation */


/*******************************************/
/****Gaussian ASAM Functions****************/
void GaussianSAMlearn();
double GaussianSAM(double in);		/* Gaussian SAM */
void GaussianSAMinit();
double GaussianSAMapprox();


/*******************************************/
/****Laplacian ASAM Functions****************/
void LaplaceSAMlearn();
double LaplaceSAM(double in);		/* Gaussian SAM */
void LaplaceSAMinit();
double LaplaceSAMapprox();



/*******************************************/
/****Triangle ASAM Functions****************/
void TriangleSAMlearn();
double TriangleSAM(double in);		/* Gaussian SAM */
void TriangleSAMinit();
double TriangleSAMapprox();



/*******************************************/
/****Cauchy ASAM Functions****************/
void CauchySAMlearn();
double CauchySAM(double in);		/* Gaussian SAM */
void CauchySAMinit();
double CauchySAMapprox();


/****Sinc ASAM Functions****************/
void SincSAMlearn();
double SincSAM(double in);		
void SincSAMinit();
double SincSAMapprox();


/*******************************************/
/****Tanh ASAM Functions****************/
void TanhSAMlearn();
double TanhSAM(double in);		/* Gaussian SAM */
void TanhSAMinit();
double TanhSAMapprox();

/*******************************************/
/******Aggregate ASAM Calls**********************/
void ASAMsInitialize();
void ASAMsLearn();
vector<double> ASAMsApprox();



#endif
