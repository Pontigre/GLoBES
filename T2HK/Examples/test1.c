#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char AEDLFILE[]="T2HK_orig.glb";
char AEDLFILE2[]="nuPRISM.glb";
char AEDLFILE3[]="Reactor2.glb";
char MYFILE[]="test1.dat";
FILE *outfile = NULL;
FILE *stream;

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
  outfile = fopen(MYFILE,"w");
  stream = stdout;

  glbInit(argv[0]);
  glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 

  theta12 = asin(sqrt(0.307));
  theta13 = asin(sqrt(0.0241));
  theta23 = 0.5;
  deltacp = 0.0;
  dm21 = 7.6e-5;
  dm32 = 2.4e-3;
  dm31 = dm32 + dm21;

  glb_params true_values = glbAllocParams();
  glb_params fit_values = glbAllocParams();
  glb_params central_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_params minimum = glbAllocParams();

  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm32);
  glbSetDensityParams(true_values,1.0,GLB_ALL);

  glbSetOscillationParameters(true_values);
  glbSetRates();

  int i;
  fprintf(outfile,"\n Oscillation probabilities in vacuum: ");
  for(i=1;i<4;i++){
    fprintf(outfile,"1->%i: %g ",i,glbVacuumProbability(1,i,+1,50,3000));}
  fprintf(outfile,"\n Oscillation probabilities in matter: ");
  for(i=1;i<4;i++){
    fprintf(outfile,"1->%i: %g ",i,glbProfileProbability(0,1,i,+1,50));}
  fprintf(outfile,"\n\n");  

  glbCopyParams(true_values,fit_values);
  glbSetOscParams(fit_values,asin(sqrt(0.03))/2,GLB_THETA_13); 

  double chi2,chi2b,chi2sum;
  chi2 = glbChiSys(fit_values,GLB_ALL,GLB_ALL);
  fprintf(outfile,"chi2 with systematics only: %g\n\n",chi2);

  chi2 = glbChiSys(fit_values,0,0);
  fprintf(outfile,"This we would have from the CP-even appearance channel only: %g\n\n",chi2);

  chi2 = glbChiSys(fit_values,GLB_ALL,0)+
         glbChiSys(fit_values,GLB_ALL,1)+
         glbChiSys(fit_values,GLB_ALL,2)+
         glbChiSys(fit_values,GLB_ALL,3);
  fprintf(outfile,"The sum over all rules gives again: %g\n\n",chi2);

  glbDefineParams(input_errors,theta12*0.1,0,0,0,dm21*0.1,0);
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetCentralValues(true_values);
  glbSetInputErrors(input_errors);

  chi2 = glbChiTheta13(fit_values,minimum,GLB_ALL); 
  fprintf(outfile,"chi2 with correlations: %g \n",chi2);
  fprintf(outfile,"Position of minimum: theta12,theta13,theta23,delta,dm21,dm32,rho\n");
  glbPrintParams(outfile,minimum);  
  fprintf(outfile,"Note that s22theta13 is unchanged/kept fixed: %g! \n\n",
	  pow(sin(2*glbGetOscParams(minimum,GLB_THETA_13)),2));

  chi2 = glbChiTheta13Delta(fit_values,minimum,GLB_ALL);
  fprintf(outfile,"chi2 with correlations other than with deltacp: %g \n\n",chi2);

  glb_projection myprojection = glbAllocProjection();
  glbDefineProjection(myprojection,GLB_FIXED,GLB_FIXED,GLB_FIXED,GLB_FREE,
   GLB_FIXED,GLB_FIXED);
  glbSetDensityProjectionFlag(myprojection,GLB_FREE,GLB_ALL);
  glbSetProjection(myprojection);
  chi2 = glbChiNP(fit_values,minimum,GLB_ALL);
  fprintf(outfile,"chi2 with correlation only with deltacp: %g \n\n",chi2);
  glbFreeProjection(myprojection);

  glbSwitchSystematics(GLB_ALL,GLB_ALL,GLB_OFF);   
  chi2 = glbChiSys(fit_values,GLB_ALL,GLB_ALL);
  glbSwitchSystematics(GLB_ALL,GLB_ALL,GLB_ON);   
  fprintf(outfile,"chi2 with statistics only: %g\n\n",chi2);

  glbDefineParams(input_errors,theta12*0.1,0,0,0,dm21*0.1,dm32/3);
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbDefineParams(central_values,theta12,theta13,theta23,deltacp,dm21,-dm32);
  glbSetDensityParams(central_values,1.0,GLB_ALL);
  glbSetCentralValues(central_values);
  glbSetInputErrors(input_errors);
  chi2=glbChiAll(central_values,minimum,GLB_ALL); 
  fprintf(outfile,"chi2 at minimum: %g \n",chi2);
  fprintf(outfile,"Position of minimum: theta12,theta13,theta23,delta,dm21,dm32,rho\n");
  glbPrintParams(outfile,minimum);  

  glbFreeParams(true_values);
  glbFreeParams(fit_values); 
  glbFreeParams(central_values);
  glbFreeParams(input_errors);
  glbFreeParams(minimum);

  return 0;
}
