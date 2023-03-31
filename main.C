// file: main.C for kbsel

// simulates evolution of 3 populations with gene flow and directional or balancing selection

// designed for ABC inference of selection on a SV locus in Timema knulli

// Time-stamp: <Tuesday, 07 December 2021, 19:54 MST -- zgompert>

#include <time.h>
#include <getopt.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "kbsel.H"

using namespace std;

gsl_rng * r;  /* global state variable for random number generator */

/* ----------------- */
/* Beginning of main */
/* ----------------- */

int main(int argc, char *argv[]) {

  time_t start = time(NULL);
  time_t end;
  int rng_seed = 0;
  int ch = 0;
  int nsims = 100;
  int i;

  string outFile = "out_kbs_abc.txt";

  // create structures for priors and parameters
  priors prior;
  params param;

  gsl_vector * pp;
  pp = gsl_vector_calloc(3); // 3 final allele frquencies

  gsl_vector_int * N; // diploid gene copy sample size
  N = gsl_vector_int_calloc(3);
  gsl_vector_int_set(N, 0, 2 * NBCTURNC);
  gsl_vector_int_set(N, 1, 2 * NBCERW);
  gsl_vector_int_set(N, 2, 2 * NBCEC);

  
  defaultPriors(&prior); // assign default priors, for now no command line argument to change defaults
  // seems easier to just recompile as these will be consistent for all simulations

  while ((ch = getopt(argc, argv, "n:o:a:b:c:A:B:C:")) != -1){
    switch(ch){
    case 'n':
      nsims = atoi(optarg);
      break;
    case 'o':
      outFile = optarg;
      break;
    case 'a':
      prior.NeAlb = atof(optarg);
      break;
    case 'b':
      prior.NeBlb = atof(optarg);
      break;
    case 'c':
      prior.NeClb = atof(optarg);
      break;
    case 'A':
      prior.NeAub = atof(optarg);
      break;
    case 'B':
      prior.NeBub = atof(optarg);
      break;
    case 'C':
      prior.NeCub = atof(optarg);
      break;
    case '?':
    default:
      cerr << "Command line failure" << endl;//usage(argv[0]);
    }
  }

  // set up gsl random number generation 
  gsl_rng_env_setup();
  r = gsl_rng_alloc (gsl_rng_default);
  srand(time(NULL));
  rng_seed = rand();
  gsl_rng_set(r, rng_seed); /* seed gsl_rng with output of rand, which
                               was seeded with result of time(NULL) */

  // open outfile
  FILE * OUT;
  OUT = fopen(outFile.c_str(), "w");
  
  // simulations
  for(i=0; i<nsims; i++){
    // sample parameter values from the relevant prior distributions
    sampleParams(&prior, &param);
    // simulate evolution via drift, selection and gene flow based
    // on the sampled parameter values
    simEvo(&param, pp);
    // summarize and write the parameter values and final allele frequencies
    writeResults(&param, pp, N, OUT);
  }

  // free memory
  gsl_vector_free(pp);

  //close outfile
  fclose(OUT);
  
  // prints run time
  end = time(NULL);
  cout << "Runtime: " << (end-start)/3600 << " hr " << (end-start)%3600/60 << " min ";
  cout << (end-start)%60 << " sec" << endl;
  return 0;
}
