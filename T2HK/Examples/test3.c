#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char AEDLFILE[]="T2HK_orig.glb";
char AEDLFILE2[]="nuPRISM.glb";
char AEDLFILE3[]="Reactor2.glb";
char MYFILE[]="test3.dat";
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

  theta12 = asin(sqrt(0.307));
  theta13 = asin(sqrt(0.0241));
  theta23 = 0.5;
  deltacp = 0.0;
  dm21 = 7.6e-5;
  dm32 = 2.4e-3;
  dm31 = dm32 + dm21;

  glb_params true_values = glbAllocParams();
  glb_params central_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_params deg_pos = glbAllocParams();

  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
  glbSetOscillationParameters(true_values);
  glbSetRates();

  glbDefineParams(central_values,theta12,theta13,theta23,deltacp,dm21,-dm31);  
  glbSetDensityParams(central_values,1.0,GLB_ALL);
  glbDefineParams(input_errors,theta12*0.1,0,0,0,dm21*0.1,dm31/3);  
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetCentralValues(central_values);
  glbSetInputErrors(input_errors);

  double CL=glbChiAll(central_values,deg_pos,GLB_ALL);
   
  printf("Position of degeneracy: s22th13=%g, deltacp=%g; Confidence level: %g \n",
    glbGetOscParams(deg_pos,GLB_THETA_13),glbGetOscParams(deg_pos,GLB_DELTA_CP),CL);

  if(CL<9.0)
  {
    double thetheta13,x,y,res;    
    for(x=-4.0;x<-2.0+0.01;x=x+2.0/50)
    for(y=0.0;y<200.0+0.01;y=y+200.0/50)
    {
        thetheta13=asin(sqrt(pow(10,x)))/2;
        glbSetOscParams(deg_pos,thetheta13,GLB_THETA_13);
        glbSetOscParams(deg_pos,y*M_PI/180.0,GLB_DELTA_CP);
    
        res=glbChiSys(deg_pos,GLB_ALL,GLB_ALL);

        fprint(outfile, "%g %g %g \n",x,y,res);
    }
  }
   
  glbFreeParams(true_values);
  glbFreeParams(central_values); 
  glbFreeParams(input_errors); 
  glbFreeParams(deg_pos); 

 return 0;
}
