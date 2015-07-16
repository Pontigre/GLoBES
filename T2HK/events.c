#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char AEDLFILE[]="T2HK.glb";
char AEDLFILE2[]="nuPRISM.glb";
char AEDLFILE3[]="Reactor2.glb";
char MYFILE[]="rates_e.dat";
char MYFILE2[]="rates_mu.dat";
FILE *outfile = NULL;
FILE *outfile2 = NULL;

/* Global parameters */
double theta12;
double theta13;
double theta23;
double deltacp;
double dm21;
double dm32;
double dm31;
glb_params true_values;

int main(int argc, char *argv[])
{ 
  glbInit(argv[0]);

  glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 
  outfile = fopen(MYFILE, "w");
  outfile2 = fopen(MYFILE2, "w");

  theta12 = asin(sqrt(0.307));
  theta13 = asin(sqrt(0.0241));
  theta23 = 0.5;
  deltacp = 0.0;
  dm21 = 7.6e-5;
  dm32 = 2.4e-3;
  dm31 = dm32 + dm21;
  int i;
  int n_bins = glbGetNumberOfBins(0);

  true_values = glbAllocParams();
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
  
  glbSetOscillationParameters(true_values);
  glbSetRates();
  double *true_rates_e = glbGetRuleRatePtr(0,5);
  double *true_rates_mu = glbGetRuleRatePtr(0,4);
  for(i=0;i<n_bins;i++){
    fprintf(outfile,"%i %g \n", i, true_rates_e[i]);
    fprintf(outfile2,"%i %g \n", i, true_rates_mu[i]);
  }

  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);

  return 0;
}

