#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char AEDLFILE[]="T2HK.glb";
char AEDLFILE2[]="nuPRISM.glb";
char AEDLFILE3[]="Reactor2.glb";
char MYFILE[]="rates.dat";
FILE *outfile = NULL;

/* Global parameters */
double theta12;
double theta13;
double theta23;
double deltacp;
double dm21;
double dm32;
glb_params true_values;

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
  deltacp = 0.0;
  dm21 = 7.6e-5;
  dm32 = 2.4e-3;

  true_values = glbAllocParams();

  /* Define "true" oscillation parameter vector */
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm32);
  glbSetDensityParams(true_values,1.0,GLB_ALL);

  /* Compute simulated data */
  glbSetOscillationParameters(true_values);
  glbSetRates();

  glbShowRuleRates(outfile, 0, 2, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG);
    
  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);

  return 0;
}

