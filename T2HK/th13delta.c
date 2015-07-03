#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char MYFILE[]="th13delta.dat";
char AEDLFILE[]="T2HK.glb";
//char AEDLFILE2[]="nuPRISM.glb";
//char AEDLFILE3[]="Reactor2.glb";
FILE *outfile = NULL;

int main(int argc, char *argv[])
{ 
  /* Initialize libglobes */
  glbInit(argv[0]);
  glbSelectMinimizer(GLB_MIN_POWELL);

  /* Initialize experiment */
  glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 
  //glbInitExperiment(AEDLFILE2,&glb_experiment_list[0],&glb_num_of_exps); 
  //glbInitExperiment(AEDLFILE3,&glb_experiment_list[0],&glb_num_of_exps); 

  /* Intitialize output */
  outfile = fopen(MYFILE, "w");
  if (outfile == NULL)
  {
    printf("Error opening output file.\n");
    return -1;
  }

  /* /\* Half the original errors *\/ */
  /* /\* NU_E_Appearance_QE *\/ */
  /* glbSetSignalErrors(GLB_ALL, 0, 5., 0.00005); */
  /* glbSetBGErrors(GLB_ALL, 0, 0.025, 0.025); */
  /* /\* NU_E_BAR_Appearance_QE *\/ */
  /* glbSetSignalErrors(GLB_ALL, 1, 5., 0.00005); */
  /* glbSetBGErrors(GLB_ALL, 1, 0.025, 0.025); */
  /* /\* NU_MU_Disapperance_QE *\/ */
  /* glbSetSignalErrors(GLB_ALL, 2, 0.00125, 0.00005); */
  /* glbSetBGErrors(GLB_ALL, 2, 0.1, 0.00005); */
  /* /\* NU_MU_BAR_Disappearance_QE *\/ */
  /* glbSetSignalErrors(GLB_ALL, 3, 0.0125, 0.00005); */
  /* glbSetBGErrors(GLB_ALL, 3, 0.1, 0.00005); */
  /* /\* NU_E_Appearance_CC *\/ */
  /* glbSetSignalErrors(GLB_ALL, 4, 0.025, 0.00005); */
  /* glbSetBGErrors(GLB_ALL, 4, 0.025, 0.00005); */
  /* /\* NU_E_BAR_Appearance_CC *\/ */
  /* glbSetSignalErrors(GLB_ALL, 5, 0.025, 0.00005); */
  /* glbSetBGErrors(GLB_ALL, 5, 0.025, 0.00005); */

  /* /\* Double the original errors *\/ */
  /* /\* NU_E_Appearance_QE *\/ */
  /* glbSetSignalErrors(GLB_ALL, 0, 20., 0.0002); */
  /* glbSetBGErrors(GLB_ALL, 0, 0.10, 0.10); */
  /* /\* NU_E_BAR_Appearance_QE *\/ */
  /* glbSetSignalErrors(GLB_ALL, 1, 20., 0.0002); */
  /* glbSetBGErrors(GLB_ALL, 1, 0.10, 0.10); */
  /* /\* NU_MU_Disapperance_QE *\/ */
  /* glbSetSignalErrors(GLB_ALL, 2, 0.05, 0.0002); */
  /* glbSetBGErrors(GLB_ALL, 2, 0.4, 0.0002); */
  /* /\* NU_MU_BAR_Disappearance_QE *\/ */
  /* glbSetSignalErrors(GLB_ALL, 3, 0.05, 0.0002); */
  /* glbSetBGErrors(GLB_ALL, 3, 0.4, 0.0002); */
  /* /\* NU_E_Appearance_CC *\/ */
  /* glbSetSignalErrors(GLB_ALL, 4, 0.10, 0.0002); */
  /* glbSetBGErrors(GLB_ALL, 4, 0.10, 0.0002); */
  /* /\* NU_E_BAR_Appearance_CC *\/ */
  /* glbSetSignalErrors(GLB_ALL, 5, 0.10, 0.0002); */
  /* glbSetBGErrors(GLB_ALL, 5, 0.10, 0.0002); */

  /* Define "true" oscillation parameters */
  double theta12 = asin(sqrt(0.307));
  double theta13 = asin(sqrt(0.0241));
  double theta23 = 0.5;
  double deltacp = 0.0;
  double dm21 = 7.6e-5;
  double dm32 = 2.4e-3;

  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_projection th13delta_projection = glbAllocProjection();
  glbDefineProjection(th13delta_projection,GLB_FIXED,GLB_FIXED,GLB_FREE,
		      GLB_FIXED,GLB_FIXED,GLB_FREE);
  glbSetDensityProjectionFlag(th13delta_projection, GLB_FIXED, GLB_ALL);
  glbSetProjection(th13delta_projection);

  for(deltacp = -M_PI/2.0; deltacp < M_PI+0.1; deltacp += M_PI/2.0){
    double x = deltacp/M_PI;
    printf("Delta_CP = %g\n", x);

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
    double this_th13, this_delta;
    double th13_lower  = asin(sqrt(0.01));
    double th13_upper  = asin(sqrt(0.04));
    double th13_steps  = 15;
    double delta_lower = -M_PI;
    double delta_upper = M_PI;
    double delta_steps = 15;
    double res;

    for(this_th13=th13_lower; this_th13<=th13_upper; this_th13+=(th13_upper-th13_lower)/th13_steps)
      {
	for(this_delta=delta_lower; this_delta<=delta_upper; this_delta+=(delta_upper-delta_lower)/delta_steps)
	  {
	    /* Set vector of test=fit values */
	    glbSetOscParams(test_values, this_th13, GLB_THETA_13);
	    glbSetOscParams(test_values, this_delta, GLB_DELTA_CP);

	    double i = sin(2*this_th13)*sin(2*this_th13);
	    double j = this_delta*180.0/M_PI;

	    /* Compute chi^2 assuming the normal mass hierarchy in the fit */
	    res = glbChiNP(test_values, NULL, GLB_ALL);
	    fprintf(outfile, "%g %g %g\n", i, j, res);
	  }
	fprintf(outfile, "\n");
      }
  }

  fclose(outfile);
  
  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);
  glbFreeProjection(th13delta_projection);

  return 0;
}

