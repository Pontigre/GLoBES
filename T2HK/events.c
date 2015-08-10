#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char AEDLFILE1[]="T2HK.glb";
char AEDLFILE2[]="nuPRISM.glb";
char AEDLFILE3[]="Reactor2.glb";
char MYFILE1[]="nuPRISM_events_e.dat";
char MYFILE2[]="nuPRISM_events_mu.dat";
char MYFILE3[]="nuPRISM_flux_e.dat";
char MYFILE4[]="nuPRISM_flux_mu.dat";
char MYFILE5[]="T2HK_events_e.dat";
char MYFILE6[]="T2HK_events_mu.dat";
char MYFILE7[]="T2HK_flux_e.dat";
char MYFILE8[]="T2HK_flux_mu.dat";
FILE *outfile1 = NULL;
FILE *outfile2 = NULL;
FILE *outfile3 = NULL;
FILE *outfile4 = NULL;
FILE *outfile5 = NULL;
FILE *outfile6 = NULL;
FILE *outfile7 = NULL;
FILE *outfile8 = NULL;

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
  
  glbInitExperiment(AEDLFILE1,&glb_experiment_list[0],&glb_num_of_exps); 
  glbInitExperiment(AEDLFILE2,&glb_experiment_list[0],&glb_num_of_exps); 

  theta12 = asin(sqrt(0.307));
  theta13 = asin(sqrt(0.0241));
  theta23 = 0.5;
  deltacp = 0.0;
  dm21 = 7.6e-5;
  dm32 = 2.4e-3;
  dm31 = dm32 + dm21;

  true_values = glbAllocParams();
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
  
  glbSetOscillationParameters(true_values);
  glbSetRates();

  int i;
  int n_bins_t2hk = glbGetNumberOfBins(0);
  int n_bins_nup = glbGetNumberOfBins(1);

  double bin_size = (1.2 - 0.4)/60;
  double e_min_t2hk = 0.0;
  double e_max_t2hk = 10.0;
  double e_step_t2hk = e_max_t2hk/n_bins_t2hk;
  double e_min_nup = 0.0;
  double e_max_nup = 10.0;
  double e_step_nup = e_max_nup/n_bins_nup;
  double E;

  /*--------------------------------------------- Initial Fluxes ------------------------------------------------------------*/

  outfile3 = fopen(MYFILE3, "w");
  outfile4 = fopen(MYFILE4, "w");
  outfile7 = fopen(MYFILE7, "w");
  outfile8 = fopen(MYFILE8, "w");

  for(E = e_min_t2hk; E<e_max_t2hk; E += e_step_t2hk){
    double flux_e = glbFlux(0, 0, E, 295, 1, 1);
    double flux_mu = glbFlux(0, 0, E, 295, 2, 1);
    fprintf(outfile7,"%g %g \n", E, 10*flux_e);
    fprintf(outfile8,"%g %g \n", E, 10*flux_mu);
  }

  for(E = e_min_nup; E<e_max_nup; E += e_step_nup){
    double flux_e_nup = glbFlux(1, 0, E, 1.0, 1, 1);
    double flux_mu_nup = glbFlux(1, 0, E, 1.0, 2, 1);
    fprintf(outfile3,"%g %g \n", E, 10*flux_e_nup);
    fprintf(outfile4,"%g %g \n", E, 10*flux_mu_nup);
  }

  fclose(outfile3);
  fclose(outfile4);
  fclose(outfile7);
  fclose(outfile8);

  /*-------------------------------------------------------- Event Rates ----------------------------------------------------*/

  outfile1 = fopen(MYFILE1, "w");
  outfile2 = fopen(MYFILE2, "w");
  outfile5 = fopen(MYFILE5, "w");
  outfile6 = fopen(MYFILE6, "w");

  double *true_rates_e_nup = glbGetRuleRatePtr(1,0);
  double *true_rates_mu_nup = glbGetRuleRatePtr(1,2);
  double *true_rates_e_t2hk = glbGetRuleRatePtr(0,0);
  double *true_rates_mu_t2hk = glbGetRuleRatePtr(0,2);
  for(i=0;i<n_bins_t2hk;i++){
    fprintf(outfile5,"%g %g \n", i*bin_size+0.4, 10*true_rates_e_t2hk[i]);
    fprintf(outfile6,"%g %g \n", i*bin_size+0.4, 10*true_rates_mu_t2hk[i]);
  }
  for(i=0;i<n_bins_nup;i++){
    fprintf(outfile1,"%g %g \n", i*bin_size+0.4, 10*true_rates_e_nup[i]);
    fprintf(outfile2,"%g %g \n", i*bin_size+0.4, 10*true_rates_mu_nup[i]);
  }

  fclose(outfile1);
  fclose(outfile2);
  fclose(outfile5);
  fclose(outfile6);

  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);

  return 0;
}

