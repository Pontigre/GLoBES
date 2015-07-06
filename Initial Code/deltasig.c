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
FILE *outfile = NULL;

/* Global parameters */
double theta12;
double theta13;
double theta23;
double deltacp;
double dm21;
double dm32;
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
  glbSelectMinimizer(GLB_MIN_POWELL);

  /* Initialize experiment */
  glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 
  //glbInitExperiment(AEDLFILE2,&glb_experiment_list[0],&glb_num_of_exps); /* nuPRISM */
  //glbInitExperiment(AEDLFILE3,&glb_experiment_list[0],&glb_num_of_exps); /* Reactor */

  /* Intitialize output */
  outfile = fopen(MYFILE, "w");
  if (outfile == NULL)
  {
    printf("Error opening output file.\n");
    return -1;
  }

  /* Mess with the errors */
  //HalfErrors();
  //DoubleErrors();

  /* Define "true" oscillation parameters */
  theta12 = asin(sqrt(0.307));
  theta13 = asin(sqrt(0.0241));
  theta23 = 0.5;
  deltacp = -M_PI;
  dm21 = 7.6e-5;
  dm32 = 2.4e-3;

  true_values = glbAllocParams();
  test_values = glbAllocParams();
  input_errors = glbAllocParams();

  /* Define "true" oscillation parameter vector */
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm32);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
 
  /* Define initial guess for the fit values */ 
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm32);  
  glbSetDensityParams(test_values,1.0,GLB_ALL);

  glbDefineParams(input_errors, theta12*0.1, 0, 0, 0, dm21*0.1, 0);
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetInputErrors(input_errors);
  glbSetCentralValues(true_values);

  /* Compute simulated data */
  glbSetOscillationParameters(true_values);
  glbSetRates();

  /* Scan the th13-delta plane */
  double this_delta;
  double delta_lower = -M_PI;
  double delta_upper = M_PI;
  double delta_steps = 40;
  double delta_step_size = (delta_upper-delta_lower)/delta_steps;
  double sig;

  for(this_delta=delta_lower; this_delta<=delta_upper; this_delta+= delta_step_size)
    {
      glbSetOscParams(test_values, this_delta, GLB_DELTA_CP);
      double i = this_delta*180.0/M_PI;
      glbSetOscParams(test_values, this_delta, GLB_DELTA_CP);
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

  return 0;
}

