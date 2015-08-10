#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char AEDLFILE[]="T2HK.glb";
char AEDLFILE2[]="nuPRISM.glb";
char AEDLFILE3[]="Reactor2.glb";
char MYFILE[]="th23dm32.dat";
char MYFILE2[]="th23dm32_05xerror.dat";
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
glb_projection th23dm32_projection;

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
  outfile2 = fopen(MYFILE2, "w");

  /* Define "true" oscillation parameters */
  theta12 = asin(sqrt(0.307));
  theta13 = asin(sqrt(0.0241));
  theta23 = 0.5;
  deltacp = 0.0;
  dm21 = 7.6e-5;
  dm32 = 2.4e-3;
  dm31 = dm32 + dm21;

  double this_th23, this_dm32;
  double th23_lower  = 0.4;
  double th23_upper  = 0.6;
  double th23_steps  = 30;
  double dm32_lower = 2.2e-3;
  double dm32_upper = 2.6e-3;
  double dm32_steps = 30;
  double res;

  true_values = glbAllocParams();
  test_values = glbAllocParams();
  input_errors = glbAllocParams();
  th23dm32_projection = glbAllocProjection();
  
  glbDefineProjection(th23dm32_projection,GLB_FIXED,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FIXED,GLB_FREE);
  glbSetDensityProjectionFlag(th23dm32_projection, GLB_FIXED, GLB_ALL);
  glbSetProjection(th23dm32_projection);


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

  for(this_th23=th23_lower; this_th23<=th23_upper; this_th23+=(th23_upper-th23_lower)/th23_steps)
    {
      printf("Theta 23 = %g\n", this_th23);
      for(this_dm32=dm32_lower; this_dm32<=dm32_upper; this_dm32+=(dm32_upper-dm32_lower)/dm32_steps)
	{
	  double new_dm31 = this_dm32 + dm21;
	  glbSetOscParams(test_values, this_th23, GLB_THETA_23);
	  glbSetOscParams(test_values, new_dm31, GLB_DM_31);
	  double i = sin(this_th23)*sin(this_th23);
	  double j = this_dm32;
	  res=glbChiNP(test_values, NULL, GLB_ALL);
	  printf("DM_32 = %g\n", j);

	  fprintf(outfile, "%g %g %g\n", i, j, res);
	}
      fprintf(outfile, "\n");
    }


  fclose(outfile);

  HalfErrors();

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

  for(this_th23=th23_lower; this_th23<=th23_upper; this_th23+=(th23_upper-th23_lower)/th23_steps)
    {
      printf("Theta 23 = %g\n", this_th23);
      for(this_dm32=dm32_lower; this_dm32<=dm32_upper; this_dm32+=(dm32_upper-dm32_lower)/dm32_steps)
	{
	  double new_dm31 = this_dm32 + dm21;
	  glbSetOscParams(test_values, this_th23, GLB_THETA_23);
	  glbSetOscParams(test_values, this_dm32, GLB_DM_31);
	  double i = sin(this_th23)*sin(this_th23);
	  double j = this_dm32;
	  res=glbChiNP(test_values, NULL, GLB_ALL);

	  fprintf(outfile2, "%g %g %g\n", i, j, res);
	}
      fprintf(outfile2, "\n");
    }


  fclose(outfile2);
  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);
  glbFreeProjection(th23dm32_projection);

  return 0;
}

