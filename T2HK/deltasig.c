#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char AEDLFILE[]="T2HK.glb";
char AEDLFILE2[]="nuPRISM.glb";
char AEDLFILE3[]="Reactor2.glb";
char MYFILE[]="deltasig.dat";
char MYFILE2[]="deltasig_inv.dat";
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
glb_params test_values;
glb_params input_errors;
glb_projection th13delta_projection;

int main(int argc, char *argv[])
{ 
  /* Initialize libglobes */
  glbInit(argv[0]);

  /* Initialize experiment */
  glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 
  //glbInitExperiment(AEDLFILE2,&glb_experiment_list[0],&glb_num_of_exps); /* nuPRISM */
  //glbInitExperiment(AEDLFILE3,&glb_experiment_list[0],&glb_num_of_exps); /* Reactor */

  /* Intitialize output */
  outfile = fopen(MYFILE, "w");


  /* Define "true" oscillation parameters */
  theta12 = asin(sqrt(0.307));
  theta13 = asin(sqrt(0.0241));
  theta23 = 0.5;
  deltacp = asin(0.0);
  dm21 = 7.6e-5;
  dm32 = 2.4e-3;
  dm31 = dm32+dm21;

  true_values = glbAllocParams();
  test_values = glbAllocParams();
  input_errors = glbAllocParams();

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

  /* Scan the delta plane */
  double this_delta;
  double delta_lower = -M_PI;
  double delta_upper = M_PI;
  double delta_steps = 72;
  double delta_step_size = (delta_upper-delta_lower)/delta_steps;
  double sig;
  printf("true delta_cp = %g \n", deltacp);
  for(this_delta=delta_lower; this_delta<=delta_upper; this_delta+= delta_step_size)
    {
      double x = sin(this_delta);
      double i = this_delta*180.0/M_PI;
      glbSetOscParams(test_values, asin(x), GLB_DELTA_CP);
      double chi2=glbChiDelta(test_values, NULL, GLB_ALL);
      glbSetOscParams(test_values, deltacp, GLB_DELTA_CP);
      double chitrue=glbChiDelta(test_values, NULL, GLB_ALL);
      sig = sqrt(abs(chi2 - chitrue));
      printf("delta = %g\n", i);
      fprintf(outfile, "%g %g\n", i, sig);
    }

  fclose(outfile);

  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);

  outfile2 = fopen(MYFILE2, "w");
  dm31 = dm21-dm32;

  true_values = glbAllocParams();
  test_values = glbAllocParams();
  input_errors = glbAllocParams();

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

  for(this_delta=delta_lower; this_delta<=delta_upper; this_delta+= delta_step_size)
    {
      double x = sin(this_delta);
      double i = this_delta*180.0/M_PI;
      glbSetOscParams(test_values, asin(x), GLB_DELTA_CP);
      double chi2=glbChiDelta(test_values, NULL, GLB_ALL);
      glbSetOscParams(test_values, deltacp, GLB_DELTA_CP);
      double chitrue=glbChiDelta(test_values, NULL, GLB_ALL);
      sig = sqrt(abs(chi2 - chitrue));
      printf("delta = %g\n", i);
      fprintf(outfile2, "%g %g\n", i, sig);
    }

  fclose(outfile2);

  return 0;
}

