#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

// macro definitions
#define VERSION "0.1 -- 07 December 2021"
#define TIME 2500
#define NBCTURNC 37
#define NBCERW 24
#define NBCEC 68



extern gsl_rng * r;

using namespace std;

// structure definitions
struct priors{
  double NeAlb; // lower bound for Ne on population A = BCTURN C
  double NeBlb; // lower bound for Ne on population B = BCE RW
  double NeClb; // lower bound for Ne on population C = BCE C
  double NeAub; // upper bound for Ne on population A = BCTURN C
  double NeBub; // upper bound for Ne on population B = BCE RW
  double NeCub; // upper bound for Ne on population C = BCE C
  double NePrec; // cond. precision for Ne 

  double NmABlb; // lower bound for Nm between A and B
  double NmAClb; // lower bound for Nm between A and C
  double NmBClb; // lower bound for Nm between B and C
  double NmABub; // upper bound for Nm between A and B
  double NmACub; // upper bound for Nm between A and C
  double NmBCub; // upper bound for Nm between B and C
  double NmPrec; // cond. precision for Nm 

  double Lslb; // log(s) lower bound
  double Lsub;
};


// structure definitions
struct params{
  int model; // boolean, 0 = both directional, 1 = bs on RW, 2 = bs on C, 3 = both bs

  double NeA; // Ne for population A = BCTURN C
  double NeB; // Ne for population B = BCE RW
  double NeC; // Ne for population C = BCE C

  double mAB; // mig. rate from A to B 
  double mBA; // mig. rate from B to A
  double mAC; // mig. rate from A to C
  double mCA; // mig. rate from C to A
  double mBC; // mig. rate from B to C
  double mCB; // mig. rate from C to B

  // selection coefficients for DS model
  // w11 = 1+s; w12 = 1+hs; w22 = 1
  double sRW;
  double sC;
  double hRW;
  double hC;

  // selection coefficients for BS model
  // w11 = 1-s1; w12 = 1; w22 = 1-s2
  double s1RW;
  double s2RW;
  double s1C;
  double s2C;

  // SV allele frequency; p is frequency of the allele favored on RW
  double pA; // allele frequency in A, B, and C, respectively
  double pB;
  double pC;
  
};

// function definitions
void defaultPriors(priors * prior);
void sampleParams(priors * prior, params * param);
void simEvo(params * param, gsl_vector * pp);
double dpDs(double p, double s, double h);
double dpBs(double p, double s1, double h2);
void writeResults(params * param, gsl_vector * pp, gsl_vector_int * N, FILE * OUT);
