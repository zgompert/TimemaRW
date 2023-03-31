#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "kbsel.H"

using namespace std;

/* Set parameters for priors */
void defaultPriors(priors * prior){

  // default lower and upper bounds on Ne
  prior->NeAlb = 100;
  prior->NeBlb = 100;
  prior->NeClb = 100;

  prior->NeAub = 1000;
  prior->NeBub = 1000;
  prior->NeCub = 1000; 

  prior->NePrec = 12;
  
  // default lower and upper bounds on Nm
  prior->NmABlb = 0.58;
  prior->NmAClb = 0.61;
  prior->NmBClb = 4.9;

  prior->NmABub = 0.69;
  prior->NmACub = 0.73;
  prior->NmBCub = 5.7;

  prior->NmPrec = 20;
  
  // default log bounds on s
  prior->Lslb = log(0.001);
  prior->Lsub = log(0.9);

  }

/* Sample parameter values from priors */
void sampleParams(priors * prior, params * param){
  double x, Nm;
  double a, b;
  double s, t;
  double mod; 

  // alpha a beta for beta
  a = 0.5 * prior->NePrec;
  b = (1-0.5) * prior->NePrec;
  // sample Ne
  x = gsl_ran_beta(r, a, b);
  param->NeA = x * (prior->NeAub - prior->NeAlb) + prior->NeAlb; // convert to rescaled beta
  x = gsl_ran_beta(r, a, b);
  param->NeB = x * (prior->NeBub - prior->NeBlb) + prior->NeBlb;
  x = gsl_ran_beta(r, a, b);
  param->NeC = x * (prior->NeCub - prior->NeClb) + prior->NeClb;

  
  // sample directional migration rates
  a = 0.5 * prior->NmPrec;
  b = (1-0.5) * prior->NmPrec;
  
  
  x = gsl_ran_beta(r, a, b);
  Nm = x * (prior->NmABub - prior->NmABlb) + prior->NmABlb; // convert to rescaled beta
  param->mAB = Nm/param->NeB; // convert to m, in this case prop of B made up of A

  x = gsl_ran_beta(r, a, b);
  Nm = x * (prior->NmABub - prior->NmABlb) + prior->NmABlb;
  param->mBA = Nm/param->NeA; 

  x = gsl_ran_beta(r, a, b);
  Nm = x * (prior->NmACub - prior->NmAClb) + prior->NmAClb; 
  param->mAC = Nm/param->NeC;
  
  x = gsl_ran_beta(r, a, b);
  Nm = x * (prior->NmACub - prior->NmAClb) + prior->NmAClb;
  param->mCA = Nm/param->NeA; 

  x = gsl_ran_beta(r, a, b);
  Nm = x * (prior->NmBCub - prior->NmBClb) + prior->NmBClb; 
  param->mBC = Nm/param->NeC; 

  x = gsl_ran_beta(r, a, b);
  Nm = x * (prior->NmBCub - prior->NmBClb) + prior->NmBClb;
  param->mCB = Nm/param->NeB; 

  // sample selection model and selection coefficients
  mod = gsl_ran_flat(r, 0, 1);
  if(mod < 0.25){
    param->model = 0;// both directional
  }
  else if(mod < 0.5){// bs on RW
    param->model = 1;
  }
  else if(mod < 0.75){// bs on C
    param->model = 2;
  }
  else{
    param->model = 3;// both bs
  }

  // both directional
  if(param->model == 0){
    param->hRW = gsl_ran_flat(r, 0, 1);
    param->hC = gsl_ran_flat(r, 0, 1);
    param->sRW = exp(gsl_ran_flat(r, prior->Lslb, prior->Lsub));
    param->sC = -1 * exp(gsl_ran_flat(r, prior->Lslb, prior->Lsub));
  }
  else if(param->model == 1){ // bs on RW
    s = exp(gsl_ran_flat(r, prior->Lslb, prior->Lsub));
    t = exp(gsl_ran_flat(r, prior->Lslb, prior->Lsub));
    if(s <= t){
      param->s1RW = s;
      param->s2RW = t;
    }
    else{
      param->s1RW = t;
      param->s2RW = s;
    }
    param->hC = gsl_ran_flat(r, 0, 1);
    param->sC = -1 * exp(gsl_ran_flat(r, prior->Lslb, prior->Lsub));
  }
  else if(param->model == 2){ // bs on C
    param->hRW = gsl_ran_flat(r, 0, 1);
    param->sRW =  exp(gsl_ran_flat(r, prior->Lslb, prior->Lsub));
    s = exp(gsl_ran_flat(r, prior->Lslb, prior->Lsub));
    t = exp(gsl_ran_flat(r, prior->Lslb, prior->Lsub));
    if(s <= t){
      param->s1C = t;
      param->s2C = s;
    }
    else{
      param->s1C = s;
      param->s2C = t;
    }
  }
  else{// bs both
    s = exp(gsl_ran_flat(r, prior->Lslb, prior->Lsub));
    t = exp(gsl_ran_flat(r, prior->Lslb, prior->Lsub));
    if(s <= t){
      param->s1RW = s;
      param->s2RW = t;
    }
    else{
      param->s1RW = t;
      param->s2RW = s;
    }
    s = exp(gsl_ran_flat(r, prior->Lslb, prior->Lsub));
    t = exp(gsl_ran_flat(r, prior->Lslb, prior->Lsub));
    if(s <= t){
      param->s1C = t;
      param->s2C = s;
    }
    else{
      param->s1C = s;
      param->s2C = t;
    }
  }
}

/* Simulte evolution with drift, selection and gene flow */
void simEvo(params * param, gsl_vector * pp){
  int k, j;
  double curP[3];
  double nexP[3];
  int N[3];

  // round Ne to int
  N[0] = floor(param->NeA + 0.5);
  N[1] = floor(param->NeB + 0.5);
  N[2] = floor(param->NeC + 0.5);

  // set initial p to 0.5 for all populations
  for(j=0; j<3; j++){
    curP[j] = 0.5;
  }

  // time loop, ensure convergence
  for(k=0; k<TIME; k++){
    
    // gene flow
    nexP[0] = curP[0] * (1 - (param->mBA + param->mCA)) + curP[1] * param->mBA + curP[2] * param->mCA;
    nexP[1] = curP[1] * (1 - (param->mAB + param->mCB)) + curP[0] * param->mAB + curP[2] * param->mCB;
    nexP[2] = curP[2] * (1 - (param->mAC + param->mBC)) + curP[0] * param->mAC + curP[1] * param->mBC;
    
    // selection
    if(param->model==0){ // directional both
      nexP[0] = nexP[0] + dpDs(nexP[0], param->sC, param->hC);
      nexP[1] = nexP[1] + dpDs(nexP[1], param->sRW, param->hRW);
      nexP[2] = nexP[2] + dpDs(nexP[2], param->sC, param->hC);
    }
    else if(param->model==1){ // bs on RW
      nexP[0] = nexP[0] + dpDs(nexP[0], param->sC, param->hC);
      nexP[1] = nexP[1] + dpBs(nexP[1], param->s1RW, param->s2RW);
      nexP[2] = nexP[2] + dpDs(nexP[2], param->sC, param->hC);
    }
    else if(param->model==2){ // bs on C
      nexP[0] = nexP[0] + dpBs(nexP[0], param->s1C, param->s2C);
      nexP[1] = nexP[1] + dpDs(nexP[1], param->sRW, param->hRW);
      nexP[2] = nexP[2] + dpBs(nexP[2], param->s1C, param->s2C);
    }
    else{ // all bs
      nexP[0] = nexP[0] + dpBs(nexP[0], param->s1C, param->s2C);
      nexP[1] = nexP[1] + dpBs(nexP[1], param->s1RW, param->s2RW);
      nexP[2] = nexP[2] + dpBs(nexP[2], param->s1C, param->s2C);
    }

    // drift
    for(j=0; j<3; j++){ // drift depends on 2Ne b/c diploids
      curP[j] = (double) gsl_ran_binomial(r, nexP[j], 2*N[j])/(2*N[j]);// / (double) N[j];
      gsl_vector_set(pp, j, curP[j]);
    }
    // cout << nexP[0] << " " << nexP[1] << " " << nexP[2] << " :: " << param->model << endl;
  }
  //cout << nexP[0] << " " << nexP[1] << " " << nexP[2] << " :: " << param->model << endl;
}

/* change by directional selection */
double dpDs(double p, double s, double h){

  double dp;

  dp = s * p * (1-p) * (p + h * (1 - 2 * p));
  return dp;

}

/* change by balancing selection */
double dpBs(double p, double s1, double s2){

  double dp;

  dp = p * (1-p) * (s2 - p * (s1 + s2));
  return dp;
}

/* write results for ABC */
void writeResults(params * param, gsl_vector * pp, gsl_vector_int * N, FILE * OUT){

  int j;
  int nx;
  double px;
    
  // print the model number
  fprintf(OUT,"%d", param->model);

  // print the migration rates
  fprintf(OUT," %.6f %.6f %.6f %.6f %.6f %.6f", param->mAB, param->mBA,
	  param->mAC, param->mCA, param->mBC, param->mCB);

  // print selection coefficients
  if(param->model==0){ // ds all
    fprintf(OUT," %.3f %.3f %.3f %.3f NA NA NA NA", param->sC, param->hC,
	    param->sRW, param->hRW);
  }
  else if (param->model==1){ // bs on RW
      fprintf(OUT," %.3f %.3f NA NA NA NA %.3f %.3f", param->sC, param->hC,
	    param->s1RW, param->s2RW);
  }
  else if (param->model==2){ // bs on C
    fprintf(OUT," NA NA %.3f %.3f %.3f %.3f NA NA", param->sRW, param->hRW,
	    param->s1C, param->s2C);
  }
  else{ // bs on both
      fprintf(OUT," NA NA NA NA %.3f %.3f %.3f %.3f", param->s1C, param->s2C,
	    param->s1RW, param->s2RW);
  }
  
  // print the allele frequency data
  for(j=0; j<3; j++){
    nx = gsl_ran_binomial(r, gsl_vector_get(pp, j), gsl_vector_int_get(N, j));
    px = (double) nx/gsl_vector_int_get(N,j);  
    fprintf(OUT," %.4f",px);
  }
 
  fprintf(OUT,"\n");
}	  
