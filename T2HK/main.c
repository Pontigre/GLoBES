#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char AEDLFILE[]="T2HK.glb";
char AEDLFILE2[]="nuPRISM.glb";
char AEDLFILE3[]="Reactor2.glb";
char MYFILE1[]="nuPRISM_init_flux_e.dat";
char MYFILE2[]="nuPRISM_init_flux_mu.dat";
char MYFILE3[]="nuPRISM_events_e.dat";
char MYFILE4[]="nuPRISM_events_mu.dat";
char MYFILE5[]="nuPRISM_deltasig.dat";
char MYFILE6[]="nuPRISM_th13delta.dat";
char MYFILE7[]="nuPRISM_th13delta_05xerror.dat";
FILE *outfile1 = NULL;
FILE *outfile2 = NULL;
FILE *outfile3 = NULL;
FILE *outfile4 = NULL;
FILE *outfile5 = NULL;
FILE *outfile6 = NULL;
FILE *outfile7 = NULL;

/* Global parameters */
double theta12;
double theta13;
double theta23;
double deltacp;
double dm21;
double dm32;
double dm31;
glb_params true_values;
glb_params test_values;
glb_params input_errors;
glb_projection th13delta_projection;

/* Set the errors at 50% the original uncertainty */
void HalfErrors()
{  
  /* NU_E_Appearance_QE */
  glbSetSignalErrors(GLB_ALL, 0, 5., 0.00005);
  glbSetBGErrors(GLB_ALL, 0, 0.025, 0.025);
  /* NU_E_BAR_Appearance_QE */
  glbSetSignalErrors(GLB_ALL, 1, 5., 0.00005);
  glbSetBGErrors(GLB_ALL, 1, 0.025, 0.025);
  /* NU_MU_Disapperance_QE */
  glbSetSignalErrors(GLB_ALL, 2, 0.00125, 0.00005);
  glbSetBGErrors(GLB_ALL, 2, 0.1, 0.00005);
  /* NU_MU_BAR_Disappearance_QE */
  glbSetSignalErrors(GLB_ALL, 3, 0.0125, 0.00005);
  glbSetBGErrors(GLB_ALL, 3, 0.1, 0.00005);
  /* NU_E_Appearance_CC */
  glbSetSignalErrors(GLB_ALL, 4, 0.025, 0.00005);
  glbSetBGErrors(GLB_ALL, 4, 0.025, 0.00005);
  /* NU_E_BAR_Appearance_CC */
  glbSetSignalErrors(GLB_ALL, 5, 0.025, 0.00005);
  glbSetBGErrors(GLB_ALL, 5, 0.025, 0.00005);
  printf("Errors halved \n");
}

void DoubleErrors()
{  
  /* NU_E_Appearance_QE */
  glbSetSignalErrors(GLB_ALL, 0, 20., 0.0002);
  glbSetBGErrors(GLB_ALL, 0, 0.10, 0.10);
  /* NU_E_BAR_Appearance_QE */
  glbSetSignalErrors(GLB_ALL, 1, 20., 0.0002);
  glbSetBGErrors(GLB_ALL, 1, 0.10, 0.10);
  /* NU_MU_Disapperance_QE */
  glbSetSignalErrors(GLB_ALL, 2, 0.05, 0.0002);
  glbSetBGErrors(GLB_ALL, 2, 0.4, 0.0002);
  /* NU_MU_BAR_Disappearance_QE */
  glbSetSignalErrors(GLB_ALL, 3, 0.05, 0.0002);
  glbSetBGErrors(GLB_ALL, 3, 0.4, 0.0002);
  /* NU_E_Appearance_CC */
  glbSetSignalErrors(GLB_ALL, 4, 0.10, 0.0002);
  glbSetBGErrors(GLB_ALL, 4, 0.10, 0.0002);
  /* NU_E_BAR_Appearance_CC */
  glbSetSignalErrors(GLB_ALL, 5, 0.10, 0.0002);
  glbSetBGErrors(GLB_ALL, 5, 0.10, 0.0002);
  printf("Errors doubled \n");
}

int main(int argc, char *argv[])
{ 
  /* Initialize libglobes */
  glbInit(argv[0]);

  /* Initialize experiment */
  //glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 
  glbInitExperiment(AEDLFILE2,&glb_experiment_list[0],&glb_num_of_exps); /* nuPRISM */
  //glbInitExperiment(AEDLFILE3,&glb_experiment_list[0],&glb_num_of_exps); /* Reactor */

  /* Intitialize output */
  outfile1 = fopen(MYFILE1, "w");
  outfile2 = fopen(MYFILE2, "w");
  outfile3 = fopen(MYFILE3, "w");
  outfile4 = fopen(MYFILE4, "w");
  outfile5 = fopen(MYFILE5, "w");
  outfile6 = fopen(MYFILE6, "w");
  outfile7 = fopen(MYFILE7, "w");

  /* Define "true" oscillation parameters */
  theta12 = asin(sqrt(0.307));
  theta13 = asin(sqrt(0.0241));
  theta23 = 0.5;
  deltacp = 0.0;
  dm21 = 7.6e-5;
  dm32 = 2.4e-3;
  dm31 = dm32 + dm21;

  /* Needed for scans */
  double e_min = 0.0;
  double e_max = 3.0;
  double dist = 295.0;
  double E = 0.0;
  int bins = glbGetNumberOfBins(0);
  int polar = +1;
  int n_bins = 100;
  int bin;
  double e_step = (e_max-e_min)/(n_bins);
  double this_th13, this_delta;
  double th13_lower  = asin(sqrt(0.01));
  double th13_upper  = asin(sqrt(0.04));
  double th13_steps  = 15;
  double delta_lower = -M_PI;
  double delta_upper = M_PI;
  double delta_steps = 15;
  double res, sig, flux;


  true_values = glbAllocParams();
  test_values = glbAllocParams();
  input_errors = glbAllocParams();
  th13delta_projection = glbAllocProjection();

  glbDefineProjection(th13delta_projection,GLB_FIXED,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FIXED,GLB_FREE);
  glbSetDensityProjectionFlag(th13delta_projection, GLB_FIXED, GLB_ALL);
  glbSetProjection(th13delta_projection);

  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
  glbSetDensityParams(true_values,1.0,GLB_ALL);

  glbSetOscillationParameters(true_values);
  glbSetRates();

  /*------------------------------------- Initial Fluxes ---------------------------------------------------*/

  printf("Initial fluxes \n");

  for(E = e_min; E<e_max; E += e_step){
    flux = glbFlux(0, 0, E, dist, 1, polar);
    fprintf(outfile1, "%g %g \n", E, flux);
  }

  for(E = e_min; E<e_max; E += e_step){
    flux = glbFlux(0, 0, E, dist, 2, polar);
    fprintf(outfile2, "%g %g \n", E, flux);
  }

  fclose(outfile1);
  fclose(outfile2);

  /*---------------------------------------- Event Rates ---------------------------------------------------*/

  printf("Event Rates \n");

  // rule 0 = NU_E_Appearance_QE
  // rule 2 = NU_MU_Disappearance_QE
  // rule 4 = NU_MU_Disappearance_CC
  // rule 5 = NU_E_Appearance_CC

  double *rates_e = glbGetRuleRatePtr(0,5);
  double *rates_mu = glbGetRuleRatePtr(0,4);

  for(bin=0;bin<bins;bin++){
    fprintf(outfile3, "%i %g \n", bin, rates_e[bin]);
    fprintf(outfile4, "%i %g \n", bin, rates_mu[bin]);
  }

  fclose(outfile3);
  fclose(outfile4);

  /*------------------------------------- Delta CP Significance ----------------------------------------------*/

  printf("Delta_cp significance \n");
  delta_steps = 40;
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm31);  
  glbSetDensityParams(test_values,1.0,GLB_ALL);
  glbDefineParams(input_errors, theta12*0.1, 0, 0, 0, dm21*0.1, 0);
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetInputErrors(input_errors);
  glbSetCentralValues(true_values);

  for(this_delta=delta_lower; this_delta<=delta_upper; this_delta+= (delta_upper-delta_lower)/delta_steps)
    {
      double x = sin(this_delta);
      double i = this_delta*180.0/M_PI;
      glbSetOscParams(test_values, asin(x), GLB_DELTA_CP);
      double chi2=glbChiDelta(test_values, NULL, GLB_ALL);
      glbSetOscParams(test_values, deltacp, GLB_DELTA_CP);
      double chitrue=glbChiDelta(test_values, NULL, GLB_ALL);
      sig = sqrt(abs(chi2 - chitrue));
      printf("delta = %g \n", i);
      fprintf(outfile5, "%g %g\n", i, sig);
    }

  fclose(outfile5);


  /*------------------------------------- theta13 - delta_cp scan ---------------------------------------------------*/

  delta_steps = 15;
  printf("Initial th13-deltacp scan \n");

  for(deltacp = -M_PI/2.0; deltacp < M_PI+0.1; deltacp += M_PI/2.0){
    double x = deltacp/M_PI;
    printf("Delta_CP = %g\n", x);

    /* Define "true" oscillation parameter vector */
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
 
    /* Define initial guess for the fit values */ 
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm31);  
    glbSetDensityParams(test_values,1.0,GLB_ALL);

    glbDefineParams(input_errors, theta12*0.1, 0, 0, 0, dm21*0.1, 0);
    glbSetDensityParams(input_errors,0.05,GLB_ALL);
    glbSetInputErrors(input_errors);
    glbSetCentralValues(true_values);

    /* Compute simulated data */
    glbSetOscillationParameters(true_values);
    glbSetRates();

    for(this_th13=th13_lower; this_th13<=th13_upper; this_th13+=(th13_upper-th13_lower)/th13_steps)
      {
	for(this_delta=delta_lower; this_delta<=delta_upper; this_delta+=(delta_upper-delta_lower)/delta_steps)
	  {
	    glbSetOscParams(test_values, this_th13, GLB_THETA_13);
	    glbSetOscParams(test_values, this_delta, GLB_DELTA_CP);
	    double i = sin(2*this_th13)*sin(2*this_th13);
	    double j = this_delta*180.0/M_PI;
	    res=glbChiNP(test_values, NULL, GLB_ALL);
	    fprintf(outfile6, "%g %g %g\n", i, j, res);
	  }
	fprintf(outfile6, "\n");
      }
  }

  fclose(outfile6);

  /*------------------------------------- theta13 - delta_cp scan - Half errors -------------------------------------*/

  printf("th13-deltacp scan with half the errors \n");

  HalfErrors();
  delta_steps = 60;
  th13_steps = 60;

  for(deltacp = -M_PI/2.0; deltacp < M_PI+0.1; deltacp += M_PI/2.0){
    double x = deltacp/M_PI;
    printf("Delta_CP = %g\n", x);

    /* Define "true" oscillation parameter vector */
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
 
    /* Define initial guess for the fit values */ 
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm31);  
    glbSetDensityParams(test_values,1.0,GLB_ALL);

    glbDefineParams(input_errors, theta12*0.1, 0, 0, 0, dm21*0.1, 0);
    glbSetDensityParams(input_errors,0.05,GLB_ALL);
    glbSetInputErrors(input_errors);
    glbSetCentralValues(true_values);

    /* Compute simulated data */
    glbSetOscillationParameters(true_values);
    glbSetRates();

    for(this_th13=th13_lower; this_th13<=th13_upper; this_th13+=(th13_upper-th13_lower)/th13_steps)
      {
	for(this_delta=delta_lower; this_delta<=delta_upper; this_delta+=(delta_upper-delta_lower)/delta_steps)
	  {
	    glbSetOscParams(test_values, this_th13, GLB_THETA_13);
	    glbSetOscParams(test_values, this_delta, GLB_DELTA_CP);
	    double i = sin(2*this_th13)*sin(2*this_th13);
	    double j = this_delta*180.0/M_PI;
	    double res=glbChiNP(test_values, NULL, GLB_ALL);
	    fprintf(outfile7, "%g %g %g\n", i, j, res);
	  }
	fprintf(outfile7, "\n");
      }
  }
  fclose(outfile7);

  glbFreeParams(test_values);
  glbFreeParams(input_errors);
  glbFreeProjection(th13delta_projection);
  glbFreeParams(true_values);
  return 0;
}

