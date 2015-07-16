#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char AEDLFILE[]="T2HK_orig.glb";
char AEDLFILE2[]="nuPRISM.glb";
char AEDLFILE3[]="Reactor2.glb";
char MYFILE[]="test2.dat";
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

int main(int argc, char *argv[])
{ 
  glbInit(argv[0]); 
  glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 

  outfile = fopen(MYFILE, "w");

  theta12 = asin(sqrt(0.8)/2);
  theta13 = asin(sqrt(0.001)/2);
  theta23 = M_PI/4;
  deltacp = M_PI/2;
  dm21 = 7.6e-5;
  dm32 = 2.4e-3;
  dm31 = dm32 + dm21;

  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_projection theta13_projection = glbAllocProjection();  

  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm31);  
  glbSetDensityParams(test_values,1.0,GLB_ALL);

  glbSetOscillationParameters(true_values);
  glbSetRates();

  glbDefineParams(input_errors,theta12*0.1,0,0,0,dm21*0.1,0);  
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetCentralValues(true_values);
  glbSetInputErrors(input_errors);


  glbDefineProjection(theta13_projection,GLB_FIXED,GLB_FIXED,GLB_FIXED,
		      GLB_FREE,GLB_FIXED,GLB_FIXED);
  glbSetDensityProjectionFlag(theta13_projection, GLB_FIXED, GLB_ALL);
  glbSetProjection(theta13_projection); 

  double thetheta13,x,res1,res2;    
  for(x=-4;x<-2+0.001;x+=2.0/50)
    {
      /* Set vector of test=fit values */
      thetheta13=asin(sqrt(pow(10,x)))/2;
      glbSetOscParams(test_values,thetheta13,GLB_THETA_13);
     
      /* Guess fit value for deltacp in order to safely find minimum */
      glbSetOscParams(test_values,200.0/2*(x+4)*M_PI/180,GLB_DELTA_CP);
 
      /* Compute Chi^2 for two-parameter correlation: minimize over deltacp only */
      res1=glbChiNP(test_values,NULL,GLB_ALL);
      
      /* Compute Chi^2 for full correlation: minimize over all but theta13 */
      res2=glbChiTheta13(test_values,NULL,GLB_ALL);
      
      fprintf(outfile, "%g %g %g \n", x, res1, res2);
    }

  glbFreeParams(true_values);
  glbFreeParams(test_values); 
  glbFreeParams(input_errors); 
  glbFreeProjection(theta13_projection);

  return 0;
}
