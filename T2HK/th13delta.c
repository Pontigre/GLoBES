#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char AEDLFILE[]="T2HK.glb";
char AEDLFILE2[]="nuPRISM.glb";
char AEDLFILE3[]="Reactor2.glb";
char MYFILE[]="th13delta.dat";
FILE *outfile = NULL;

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
  deltacp = 0.0;
  dm21 = 7.6e-5;
  dm32 = 2.4e-3;
  dm31 = dm32 + dm21;

  double this_th13, this_delta;
  double th13_lower  = asin(sqrt(0.01));
  double th13_upper  = asin(sqrt(0.04));
  double th13_steps  = 30;
  double delta_lower = -M_PI;
  double delta_upper = M_PI;
  double delta_steps = 40;
  double res;

  true_values = glbAllocParams();
  test_values = glbAllocParams();
  input_errors = glbAllocParams();
  th13delta_projection = glbAllocProjection();
  
  glbDefineProjection(th13delta_projection,GLB_FIXED,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FIXED,GLB_FREE);
  glbSetDensityProjectionFlag(th13delta_projection, GLB_FIXED, GLB_ALL);
  glbSetProjection(th13delta_projection);

  //for(deltacp = -M_PI/2.0; deltacp < M_PI+0.1; deltacp += M_PI/2.0){
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
	printf("Theta 13 = %g\n", this_th13);
	for(this_delta=delta_lower; this_delta<=delta_upper; this_delta+=(delta_upper-delta_lower)/delta_steps)
	  {
	    glbSetOscParams(test_values, this_th13, GLB_THETA_13);
	    glbSetOscParams(test_values, this_delta, GLB_DELTA_CP);
	    double i = sin(2*this_th13)*sin(2*this_th13);
	    double j = this_delta*180.0/M_PI;
	    res=glbChiNP(test_values, NULL, GLB_ALL);

	    fprintf(outfile, "%g %g %g\n", i, j, res);
	  }
	fprintf(outfile, "\n");
      }
    //}

  fclose(outfile);

  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);
  glbFreeProjection(th13delta_projection);

  return 0;
}

