#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */

/* Output file */
char MYFILE[]="th13delta.dat";
char MYFILE2[]="deltasig.dat";
char AEDLFILE[]="T2HK.glb";
char AEDLFILE2[]="nuPRISM.glb";
FILE *outfile = NULL;
FILE *outfile2 = NULL;

int main(int argc, char *argv[])
{ 
  /* Initialize libglobes */
  glbInit(argv[0]); 
  glbSelectMinimizer(GLB_MIN_POWELL);

  /* Initialize experiment */
  glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 
  glbInitExperiment(AEDLFILE2,&glb_experiment_list[0],&glb_num_of_exps); 

  /* Intitialize output */
  outfile = fopen(MYFILE, "w");
  outfile2 = fopen(MYFILE2, "w");

  /* Oscillation parameters as per Table XVIII in HK Paper */
  const double theta12 = asin(sqrt(0.8704))/2.0;
  const double theta13 = asin(sqrt(0.10))/2.0;
  const double theta23 = asin(sqrt(0.5));
  const double deltacp = 0.0;
  const double dm21 = 7.6e-5;
  const double dm32 = 2.4e-3;

  /* Define "true" oscillation parameters (cf. hep-ph/0405172v5) */
  /*double theta12 = asin(sqrt(0.307));
  double theta13 = asin(sqrt(0.0241));
  double theta23 = 0.5;
  double deltacp = asin(0.0);
  double dm21 = 7.6e-5;
  double dm32 = 2.41e-3;*/

  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_projection th13delta_projection = glbAllocProjection();

  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm32);
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm32);
  glbDefineParams(input_errors, theta12*0.1, 0, 0, 0, dm21*0.1, 0);

  glbSetDensityParams(true_values,1.0,GLB_ALL);
  glbSetDensityParams(test_values,1.0,GLB_ALL);
  glbSetDensityParams(input_errors,0.05,GLB_ALL);

  glbSetInputErrors(input_errors);
  glbSetCentralValues(true_values);

  glbDefineProjection(th13delta_projection,GLB_FIXED,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FIXED,GLB_FREE);
  glbSetDensityProjectionFlag(th13delta_projection, GLB_FIXED, GLB_ALL);
  glbSetProjection(th13delta_projection);  
  
  /* Compute simulated data */
  glbSetOscillationParameters(true_values);
  glbSetRates();

  /* ----------------------------------------- SCAN OVER DELTA CP -------------------------------------------------*/

 /* Scan over delta */
  double this_dcp = 0.0;
  double dcp_lower = -1.0*M_PI;
  double dcp_upper = 1.0*M_PI;
  double dcp_steps = 20;
  double dcp_step_size = (dcp_upper-dcp_lower)/dcp_steps;
  double chi_true = 0.0;
  double chi2 = 0.0;
  double sig = 0.0;
  
  for(this_dcp=dcp_lower; this_dcp<=dcp_upper; this_dcp+=dcp_step_size)
    { 
      double j = this_dcp/M_PI;
      printf("delta = %g\n", j);

      glbSetOscParams(test_values, deltacp, GLB_DELTA_CP);
      chi_true = glbChiDelta(true_values, NULL, GLB_ALL);

      glbSetOscParams(test_values, this_dcp, GLB_DELTA_CP);
      chi2 = glbChiDelta(test_values, NULL, GLB_ALL);

      sig = sqrt(abs(chi2-chi_true));
      fprintf(outfile2, "%g %g\n", j, sig);
    }
  
  fclose(outfile2);

  /* ------------------------------ FOR DELTA_CP = -90 TO 180 DEG - TH13-DELTA PLANE -------------------------------*/

  /* Needed to scan the th13-delta plane */
  /* double this_th13, this_delta;
  double th13_lower  = asin(sqrt(0.01))/2.0;
  double th13_upper  = asin(sqrt(0.8))/2.0;
  double th13_steps  = 20;
  double th13_step_size = (th13_upper-th13_lower)/th13_steps;
  double delta_lower = -M_PI;
  double delta_upper = M_PI;
  double delta_steps = 20;
  double delta_step_size = (delta_upper-delta_lower)/delta_steps;
  double res;

  for(this_dcp = -M_PI/2.0; this_dcp < M_PI+0.1; this_dcp += M_PI/2.0)
    {
  printf("Delta_CP = %g\n", this_dcp/M_PI);
 
      glbDefineParams(true_values,theta12,theta13,theta23,this_dcp,dm21,dm32);
      glbSetDensityParams(true_values,1.0,GLB_ALL);

      glbDefineParams(test_values,theta12,theta13,theta23,this_dcp,dm21,dm32);  
      glbSetDensityParams(test_values,1.0,GLB_ALL);

      glbSetOscillationParameters(true_values);
      glbSetRates();

      for(this_th13=th13_lower; this_th13<=th13_upper; this_th13+=th13_step_size)
	{
	  printf("th13 = %g\n", this_th13);
	  for(this_delta=delta_lower; this_delta<=delta_upper; this_delta+=delta_step_size)
	    {
	      glbSetOscParams(test_values, this_th13, GLB_THETA_13);
	      glbSetOscParams(test_values, this_delta, GLB_DELTA_CP);

	      double i = sin(2*this_th13)*sin(2*this_th13);
	      double j = this_delta*180.0/M_PI;

	      res = glbChiNP(test_values, NULL, GLB_ALL);
	      fprintf(outfile, "%g %g %g\n", i, j, res);
	    }
	  fprintf(outfile, "\n");
	}
      }
  */
  fclose(outfile);
  /* ----------------------------------------------END-------------------------------------------------------------*/

  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);
  glbFreeProjection(th13delta_projection);
  
  return 0;
}
