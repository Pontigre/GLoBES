#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */

/* Output file */
char MYFILE[]="deltasig.dat";
char MYFILE2[]="th13delta.dat";
char MYFILE3[]="deltasig_errors2x.dat";
char MYFILE4[]="th13delta_errors2x.dat";
char MYFILE5[]="deltasig_errors05x.dat";
char MYFILE6[]="th13delta_errors05x.dat";
char AEDLFILE[]="T2HK.glb";
char AEDLFILE2[]="nuPRISM.glb";
FILE *outfile = NULL;
FILE *outfile2 = NULL;
FILE *outfile3 = NULL;
FILE *outfile4 = NULL;
FILE *outfile5 = NULL;
FILE *outfile6 = NULL;

int main(int argc, char *argv[])
{ 
  int th13delta = 1;
  int deltasig = 0;
  int th13delta_errors2x = 0;
  int deltasig_errors2x = 0;
  int th13delta_errors05x = 0;
  int deltasig_errors05x = 0;

  if(th13delta){
    printf("th13delta = true \n");}
  else{
    printf("th13delta = false \n");}

  if(deltasig){
    printf("deltasig = true \n");}
  else{
    printf("deltasig = false \n");}

  if(th13delta_errors2x){
    printf("th13delta_errors2x = true \n");}
  else{
    printf("th13delta_errors2x = false \n");}

  if(deltasig_errors2x){
    printf("deltasig_errors2x = true \n");}
  else{
    printf("deltasig_errors2x = false \n");}

  if(th13delta_errors05x){
    printf("th13delta_errors2x = true \n");}
  else{
    printf("th13delta_errors2x = false \n");}

  if(deltasig_errors05x){
    printf("deltasig_errors0.5x = true \n");}
  else{
    printf("deltasig_errors0.5x = false \n");}

  /* Initialize libglobes */
  glbInit(argv[0]); 
  glbSelectMinimizer(GLB_MIN_POWELL);

  /* Initialize experiment */
  glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 
  glbInitExperiment(AEDLFILE2,&glb_experiment_list[0],&glb_num_of_exps); 

  /* Define "true" oscillation parameters (cf. hep-ph/0405172v5) */
  /* double theta12 = asin(sqrt(0.307)); */
  /* double theta13 = asin(sqrt(0.0241)); */
  /* double theta23 = 0.5; */
  /* double deltacp = asin(0.0); */
  /* double dm21 = 7.6e-5; */
  /* double dm32 = 2.41e-3; */

  /* Oscillation parameters as per Table XVIII in HK Paper */
  const double theta12 = asin(sqrt(0.8704))/2.0;
  const double theta13 = asin(sqrt(0.10))/2.0;
  const double theta23 = asin(sqrt(0.5));
  const double deltacp = 0.0;
  const double dm21 = 7.6e-5;
  const double dm32 = 2.4e-3;

  /* Parameters for delta cp scan */
  double this_dcp = 0.0;
  double dcp_lower = -1.0*M_PI;
  double dcp_upper = 1.0*M_PI;
  double dcp_steps = 20;
  double dcp_step_size = (dcp_upper-dcp_lower)/dcp_steps;
  double chi_true = 0.0;
  double chi2 = 0.0;
  double sig = 0.0;

  /* Parameters for th13-deltacp scan */
  double this_th13, this_delta;
  double th13_lower  = asin(sqrt(0.01))/2.0;
  double th13_upper  = asin(sqrt(0.8))/2.0;
  double th13_steps  = 20;
  double th13_step_size = (th13_upper-th13_lower)/th13_steps;
  double delta_lower = -M_PI;
  double delta_upper = M_PI;
  double delta_steps = 20;
  double delta_step_size = (delta_upper-delta_lower)/delta_steps;
  double res;
  
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

  /*   glbDefineProjection(glb_projection in, theta12, theta13, theta23, delta, dm21, dm32); */
  glbDefineProjection(th13delta_projection,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FIXED);
  glbSetDensityProjectionFlag(th13delta_projection, GLB_FIXED, GLB_ALL);
  glbSetProjection(th13delta_projection);  
  
  /* Compute simulated data */
  glbSetOscillationParameters(true_values);
  glbSetRates();

  /* ----------------------------------------- SCAN OVER DELTA CP -------------------------------------------------*/

  if(deltasig){

    outfile = fopen(MYFILE, "w");

    for(this_dcp=dcp_lower; this_dcp<=dcp_upper; this_dcp+=dcp_step_size)
      {
	double j = this_dcp/M_PI;
	printf("delta = %g\n", j);

	glbSetOscParams(test_values, deltacp, GLB_DELTA_CP);
	chi_true = glbChiDelta(true_values, NULL, GLB_ALL);

	glbSetOscParams(test_values, this_dcp, GLB_DELTA_CP);
	chi2 = glbChiDelta(test_values, NULL, GLB_ALL);

	sig = sqrt(abs(chi2-chi_true));
	fprintf(outfile, "%g %g\n", j, sig);
      }
  
    fclose(outfile);
  }
  /* ------------------------------ FOR DELTA_CP = -90 TO 180 DEG - TH13-DELTA PLANE -------------------------------*/
 
  if(th13delta){
    outfile2 = fopen(MYFILE2, "w");
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
		fprintf(outfile2, "%g %g %g\n", i, j, res);
	      }
	    fprintf(outfile2, "\n");
	  }
      }
    fclose(outfile2);
  }

  /* ---------------------------------------- DOUBLE THE ERRORS ----------------------------------------------------*/
 
  /* Double the errors */
  if(deltasig_errors2x || th13delta_errors2x){
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

    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm32);
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm32);
    glbDefineParams(input_errors, theta12*0.1, 0, 0, 0, dm21*0.1, 0);

    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbSetDensityParams(test_values,1.0,GLB_ALL);
    glbSetDensityParams(input_errors,0.05,GLB_ALL);

    glbSetInputErrors(input_errors);
    glbSetCentralValues(true_values);

    glbDefineProjection(th13delta_projection,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FIXED);
    glbSetDensityProjectionFlag(th13delta_projection, GLB_FIXED, GLB_ALL);
    glbSetProjection(th13delta_projection);
  
    /* Compute simulated data */
    glbSetOscillationParameters(true_values);
    glbSetRates();
  }
  /* ----------------------------------------- SCAN OVER DELTA CP -------------------------------------------------*/

  if(deltasig_errors2x){
    outfile3 = fopen(MYFILE3, "w");

    for(this_dcp=dcp_lower; this_dcp<=dcp_upper; this_dcp+=dcp_step_size)
      {
	double j = this_dcp/M_PI;
	printf("delta = %g\n", j);

	glbSetOscParams(test_values, deltacp, GLB_DELTA_CP);
	chi_true = glbChiDelta(true_values, NULL, GLB_ALL);

	glbSetOscParams(test_values, this_dcp, GLB_DELTA_CP);
	chi2 = glbChiDelta(test_values, NULL, GLB_ALL);

	sig = sqrt(abs(chi2-chi_true));
	fprintf(outfile3, "%g %g\n", j, sig);
      }
    fclose(outfile3);
  }

  /* ------------------------------ FOR DELTA_CP = -90 TO 180 DEG - TH13-DELTA PLANE -------------------------------*/

  if(th13delta_errors2x){
    outfile4 = fopen(MYFILE4, "w");

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
		fprintf(outfile4, "%g %g %g\n", i, j, res);
	      }
	    fprintf(outfile4, "\n");
	  }
      }
    fclose(outfile4);
  }


  /* ---------------------------------------- HALVE THE ERRORS ----------------------------------------------------*/
 
  /* Halve the errors */
  if(deltasig_errors05x || th13delta_errors05x){
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

    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm32);
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm32);
    glbDefineParams(input_errors, theta12*0.1, 0, 0, 0, dm21*0.1, 0);

    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbSetDensityParams(test_values,1.0,GLB_ALL);
    glbSetDensityParams(input_errors,0.05,GLB_ALL);

    glbSetInputErrors(input_errors);
    glbSetCentralValues(true_values);

    glbDefineProjection(th13delta_projection,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FIXED);
    glbSetDensityProjectionFlag(th13delta_projection, GLB_FIXED, GLB_ALL);
    glbSetProjection(th13delta_projection);
  
    /* Compute simulated data */
    glbSetOscillationParameters(true_values);
    glbSetRates();
  }

  if(deltasig_errors05x){
    outfile5 = fopen(MYFILE5, "w");

    for(this_dcp=dcp_lower; this_dcp<=dcp_upper; this_dcp+=dcp_step_size)
      {
	double j = this_dcp/M_PI;
	printf("delta = %g\n", j);

	glbSetOscParams(test_values, deltacp, GLB_DELTA_CP);
	chi_true = glbChiDelta(true_values, NULL, GLB_ALL);

	glbSetOscParams(test_values, this_dcp, GLB_DELTA_CP);
	chi2 = glbChiDelta(test_values, NULL, GLB_ALL);

	sig = sqrt(abs(chi2-chi_true));
	fprintf(outfile5, "%g %g\n", j, sig);
      }
    fclose(outfile5);
  }

  if(th13delta_errors05x){
    outfile6 = fopen(MYFILE6, "w");

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
		fprintf(outfile6, "%g %g %g\n", i, j, res);
	      }
	    fprintf(outfile6, "\n");
	  }
      }
    fclose(outfile6);
  }

  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);
  glbFreeProjection(th13delta_projection);
  
  return 0;
}
