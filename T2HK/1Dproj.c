#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char AEDLFILE[]="T2HK.glb";
char AEDLFILE2[]="nuPRISM.glb";
char AEDLFILE3[]="Reactor2.glb";
char MYFILE1[]="th13.dat";
char MYFILE2[]="th23.dat";
FILE *outfile1 = NULL;
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
  glbInit(argv[0]);
  glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 
  glbInitExperiment(AEDLFILE2,&glb_experiment_list[0],&glb_num_of_exps); 

  outfile1 = fopen(MYFILE1, "w");
  outfile2 = fopen(MYFILE2, "w");

  theta12 = asin(sqrt(0.307));
  theta13 = asin(sqrt(0.0241));
  theta23 = 0.5;
  deltacp = 0.0;
  dm21 = 7.6e-5;
  dm32 = 2.4e-3;
  dm31 = dm32 + dm21;

  true_values = glbAllocParams();
  test_values = glbAllocParams();
  input_errors = glbAllocParams();
  th13delta_projection = glbAllocProjection();

  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
 
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm31);  
  glbSetDensityParams(test_values,1.0,GLB_ALL);

  glbDefineParams(input_errors, theta12*0.1, 0, 0, 0, dm21*0.1, 0);
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetInputErrors(input_errors);
  glbSetCentralValues(true_values);

  glbDefineProjection(th13delta_projection,GLB_FIXED,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FIXED,GLB_FREE);
  glbSetDensityProjectionFlag(th13delta_projection, GLB_FIXED, GLB_ALL);
  glbSetProjection(th13delta_projection);

  glbSetOscillationParameters(true_values);
  glbSetRates();

  double this_th13, this_th23, chi2, i;
  double th23_lower = asin(sqrt(0.22));
  double th23_upper = asin(sqrt(0.24));
  double th23_steps = 40;
  double th23_step_size = (th23_upper-th23_lower)/th23_steps;
  double th13_lower  = asin(sqrt(0.01));
  double th13_upper  = asin(sqrt(0.04));
  double th13_steps  = 40;
  double th13_step_size = (th13_upper-th13_lower)/th13_steps;

  /* for(this_th13=th13_lower; this_th13<=th13_upper; this_th13+= th13_step_size) */
  /*   { */
  /*     i = sin(2*this_th13)*sin(2*this_th13); */
  /*     glbSetOscParams(test_values, this_th13, GLB_THETA_13); */
  /*     chi2=glbChiTheta13(test_values, NULL, GLB_ALL); */
  /*     printf("th13 = %g\n", i); */
  /*     fprintf(outfile1, "%g %g\n", i, chi2); */
  /*   } */

  /* fclose(outfile1); */

  for(this_th23=th23_lower; this_th23<=th23_upper; this_th23+= th23_step_size)
    {
      i = sin(this_th23)*sin(this_th23);
      glbSetOscParams(test_values, this_th23, GLB_THETA_23);
      chi2=glbChiTheta23(test_values, NULL, GLB_ALL);
      printf("th23 = %g\n", i);
      fprintf(outfile2, "%g %g\n", i, chi2);
    }

  fclose(outfile2);
  
  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);
  glbFreeProjection(th13delta_projection);

  return 0;
}

