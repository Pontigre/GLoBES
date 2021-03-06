#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */

/* Output file */
char MYFILE[]="th13delta.dat";
char MYFILE2[]="deltasig.dat";
char AEDLFILE[]="T2HK.glb";
FILE *outfile = NULL;
FILE *outfile2 = NULL;

int main(int argc, char *argv[])
{ 
  /* Initialize libglobes */
  glbInit(argv[0]); 
  glbSelectMinimizer(GLB_MIN_POWELL);

  /* Initialize experiment */
  glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 

  /* Intitialize output */
  outfile = fopen(MYFILE, "w");
  if (outfile == NULL)
    {
      printf("Error opening output file.\n");
      return -1;
    }
  outfile2 = fopen(MYFILE2, "w");
  if (outfile == NULL)
    {
      printf("Error opening output file.\n");
      return -1;
    }

  /* Define "true" oscillation parameters (cf. hep-ph/0405172v5) */
  double theta12 = asin(sqrt(0.307));
  double theta13 = asin(sqrt(0.0241));
  double theta23 = 0.5;
  double deltacp = 0.0;
  double sdm = 7.6e-5;
  double ldm = 2.41e-3;
  
  glb_params true_values = glbAllocParams();
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
 
  glb_params test_values = glbAllocParams();
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  
  glbSetDensityParams(test_values,1.0,GLB_ALL);

  glb_params input_errors = glbAllocParams();
  glbDefineParams(input_errors, theta12*0.1, 0, 0, 0, sdm*0.1, 0);
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetInputErrors(input_errors);
  glbSetCentralValues(true_values);

  glb_projection th13delta_projection = glbAllocProjection();
  glbDefineProjection(th13delta_projection,GLB_FIXED,GLB_FIXED,GLB_FREE,
		      GLB_FIXED,GLB_FIXED,GLB_FREE);
  glbSetDensityProjectionFlag(th13delta_projection, GLB_FIXED, GLB_ALL);
  glbSetProjection(th13delta_projection);  
  
  /* Compute simulated data */
  glbSetOscillationParameters(true_values);
  glbSetRates();

  /* Scan the th13-delta plane */
  double this_th13, this_delta;
  double th13_lower  = asin(sqrt(0.01))/2.0;
  double th13_upper  = asin(sqrt(0.4))/2.0;
  double th13_steps  = 20;
  double delta_lower = 0.0;
  double delta_upper = 2*M_PI;
  double delta_steps = 20;
  double res;
  double sig;

  /*
  for(this_th13=th13_lower; this_th13<=th13_upper; this_th13+=(th13_upper-th13_lower)/th13_steps)
    {
      for(this_delta=delta_lower; this_delta<=delta_upper; this_delta+=(delta_upper-delta_lower)/delta_steps)
	{
	  glbSetOscParams(test_values, this_th13, GLB_THETA_13);
	  glbSetOscParams(test_values, this_delta, GLB_DELTA_CP);

	  double i = sin(2*this_th13)*sin(2*this_th13);
	  double j = this_delta/M_PI;

	  res = glbChiNP(test_values, NULL, GLB_ALL);
	  fprintf(outfile, "%g %g %g\n", i, j, res);
	}
      fprintf(outfile, "\n");
      printf("theta13 = %g\n", this_th13);
    }
  */

  fclose(outfile);

  for(this_delta=delta_lower; this_delta<=delta_upper; this_delta+=(delta_upper-delta_lower)/delta_steps)
    { 
      double j = this_delta/M_PI;
      printf("delta = %g\n", j);

      glbSetOscParams(test_values, deltacp, GLB_DELTA_CP);
      double chi_true = glbChiSys(test_values, GLB_ALL, GLB_ALL);

      glbSetOscParams(test_values, this_delta, GLB_DELTA_CP);
      double chi2 = glbChiSys(test_values, GLB_ALL, GLB_ALL);
      
      double sig = sqrt(abs(chi2-chi_true));

      fprintf(outfile2, "%g %g\n", j, sig);
      printf("true chi2 = %g, test chi2 = %g, sig = %g\n", chi_true, chi2, sig);
    }    
  fclose(outfile2);
      
  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);
  glbFreeProjection(th13delta_projection);
  
  return 0;
}

