#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <globes/globes.h>   /* GLoBES library */

/* Files */
char AEDLFILE[]="T2HK_orig.glb";
char AEDLFILE2[]="nuPRISM.glb";
char AEDLFILE3[]="Reactor2.glb";
char MYFILE1[]="test1.dat";
char MYFILE2[]="test2.dat";
char MYFILE3[]="test3.dat";
char MYFILE4[]="test4a.dat";
char MYFILE5[]="test4b.dat";
char MYFILE6[]="test4c.dat";
char MYFILE7[]="test4d.dat";
FILE *outfile1 = NULL;
FILE *outfile2 = NULL;
FILE *outfile3 = NULL;
FILE *outfile4 = NULL;
FILE *outfile5 = NULL;
FILE *outfile6 = NULL;
FILE *outfile7 = NULL;

/* Global parameters */
double theta12;
double theta13;
double theta23;
double deltacp;
double dm21;
double dm32;
double dm31;
glb_params true_values;


double Sign(double n)
{
 if(n>=0) return +1.0;
 else return -1.0;   
}

/* Set simulated rates */
void SetOscParams(double thedm31)
{
 glb_params true_values = glbAllocParams();
 glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,thedm31);
 glbSetDensityParams(true_values,1.0,GLB_ALL);
 glbSetOscillationParameters(true_values);
 glbSetRates();
 glbFreeParams(true_values);  
}

/* Calculate chi^2 with systematics */
double CalcSystematics(double thedm31,double thex)
{
  SetOscParams(thedm31);
  double thetheta13=asin(sqrt(pow(10,thex)))/2;
  glb_params test_values = glbAllocParams();
  glbDefineParams(test_values,theta12,thetheta13,theta23,0.0,dm21,thedm31);
  glbSetDensityParams(test_values,1.0,GLB_ALL);
  double res=glbChiSys(test_values,GLB_ALL,GLB_ALL);
  glbFreeParams(test_values);  
  return res;
}

/* Calculate chi^2 with statistics only */
double CalcNoSystematics(double thedm31,double thex)
{
 int rules=glbGetNumberOfRules(0);
 int j;
 glbSwitchSystematics(GLB_ALL,GLB_ALL,GLB_OFF);
 double res=CalcSystematics(thedm31,thex);
 glbSwitchSystematics(GLB_ALL,GLB_ALL,GLB_ON);
 return res; 
}

/* Calculate chi^2 with correlations */
double CalcProjection(double thedm31,double thex,glb_params start_vector)
{
  SetOscParams(thedm31);

  double thetheta13=asin(sqrt(pow(10,thex)))/2;
  double thesign=Sign(glbGetOscParams(start_vector,GLB_DM_ATM));

  glb_params input_errors = glbAllocParams();
  glb_params central_values = glbAllocParams();
  glb_params minimum = glbAllocParams();

  glbDefineParams(input_errors,theta12*0.1,0,0,0,dm21*0.1,thedm31/3);  
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbDefineParams(central_values,theta12,theta13,theta23,deltacp,dm21,thedm31*thesign);
  glbSetDensityParams(central_values,1.0,GLB_ALL);
  glbSetCentralValues(central_values);
  glbSetInputErrors(input_errors);

  start_vector=glbSetOscParams(start_vector,thetheta13,GLB_THETA_13);
  double res=glbChiTheta13(start_vector,minimum,GLB_ALL);
  
  glbSetOscParams(start_vector,glbGetOscParams(minimum,GLB_DELTA_CP),GLB_DELTA_CP);
 
  glbFreeParams(input_errors);
  glbFreeParams(central_values);
  glbFreeParams(minimum);
  
  return res; 
}

/* Find degeneracy */
double FindDeg(glb_params deg_pos,double thedm31)
{
  SetOscParams(thedm31);

  glb_params input_errors = glbAllocParams();
  glb_params central_values = glbAllocParams();

  glbDefineParams(central_values,theta12,theta13,theta23,deltacp,dm21,-thedm31);  
  glbSetDensityParams(central_values,1.0,GLB_ALL);
  glbDefineParams(input_errors,theta12*0.1,0,0,0,dm21*0.1,thedm31/3);  
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetCentralValues(central_values);
  glbSetInputErrors(input_errors);
  double CL=glbChiAll(central_values,deg_pos,GLB_ALL); 

  glbFreeParams(input_errors);
  glbFreeParams(central_values);
  return CL;
}

int main(int argc, char *argv[])
{ 
  outfile1 = fopen(MYFILE1,"w");
  outfile2 = fopen(MYFILE2,"w");
  outfile3 = fopen(MYFILE3,"w");
  outfile4 = fopen(MYFILE4,"w");
  outfile5 = fopen(MYFILE5,"w");
  outfile6 = fopen(MYFILE6,"w");
  outfile7 = fopen(MYFILE7,"w");

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
  glb_params test_values = glbAllocParams();
  glb_params fit_values = glbAllocParams();
  glb_params central_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_params minimum = glbAllocParams();
  glb_projection myprojection = glbAllocProjection();
  glb_projection theta13_projection = glbAllocProjection();  
  glb_params deg_pos = glbAllocParams();
  glb_params original1 = glbAllocParams();
  glb_params original2 = glbAllocParams();
  glb_params degeneracy1 = glbAllocParams();
  glb_params degeneracy2 = glbAllocParams();


  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm32);
  glbSetDensityParams(true_values,1.0,GLB_ALL);

  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm31);  
  glbSetDensityParams(test_values,1.0,GLB_ALL);

  glbDefineParams(input_errors,theta12*0.1,0,0,0,dm21*0.1,0);

  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetInputErrors(input_errors);
  glbSetCentralValues(true_values);

  glbDefineProjection(theta13_projection,GLB_FIXED,GLB_FIXED,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FIXED);
  glbSetDensityProjectionFlag(theta13_projection,GLB_FIXED,GLB_ALL);
  glbSetProjection(theta13_projection); 

  glbDefineProjection(myprojection,GLB_FIXED,GLB_FIXED,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FIXED);
  glbSetDensityProjectionFlag(myprojection,GLB_FREE,GLB_ALL);
  glbSetProjection(myprojection);

  glbSetOscillationParameters(true_values);
  glbSetRates();

  /*---------------------------- Test 1 --------------------------*/

  printf("Test 1\n");
  int i;
  fprintf(outfile1,"\n Oscillation probabilities in vacuum: ");
  for(i=1;i<4;i++){
    fprintf(outfile1,"1->%i: %g ",i,glbVacuumProbability(1,i,+1,50,3000));}
  fprintf(outfile1,"\n Oscillation probabilities in matter: ");
  for(i=1;i<4;i++){
    fprintf(outfile1,"1->%i: %g ",i,glbProfileProbability(0,1,i,+1,50));}
  fprintf(outfile1,"\n\n");  

  glbCopyParams(true_values,fit_values);
  glbSetOscParams(fit_values,asin(sqrt(0.03))/2,GLB_THETA_13); 

  double chi2,chi2b,chi2sum;
  chi2 = glbChiSys(fit_values,GLB_ALL,GLB_ALL);
  fprintf(outfile1,"chi2 with systematics only: %g\n\n",chi2);

  chi2 = glbChiSys(fit_values,0,0);
  fprintf(outfile1,"This we would have from the CP-even appearance channel only: %g\n\n",chi2);

  chi2 = glbChiSys(fit_values,GLB_ALL,0)+
         glbChiSys(fit_values,GLB_ALL,1)+
         glbChiSys(fit_values,GLB_ALL,2)+
         glbChiSys(fit_values,GLB_ALL,3);
  fprintf(outfile1,"The sum over all rules gives again: %g\n\n",chi2);

  chi2 = glbChiTheta13(fit_values,minimum,GLB_ALL); 
  fprintf(outfile1,"chi2 with correlations: %g \n",chi2);
  fprintf(outfile1,"Position of minimum: theta12,theta13,theta23,delta,dm21,dm32,rho\n");
  glbPrintParams(outfile1,minimum);  
  fprintf(outfile1,"Note that s22theta13 is unchanged/kept fixed: %g! \n\n",
	  pow(sin(2*glbGetOscParams(minimum,GLB_THETA_13)),2));

  chi2 = glbChiTheta13Delta(fit_values,minimum,GLB_ALL);
  fprintf(outfile1,"chi2 with correlations other than with deltacp: %g \n\n",chi2);


  chi2 = glbChiNP(fit_values,minimum,GLB_ALL);
  fprintf(outfile1,"chi2 with correlation only with deltacp: %g \n\n",chi2);

  glbSwitchSystematics(GLB_ALL,GLB_ALL,GLB_OFF);   
  chi2 = glbChiSys(fit_values,GLB_ALL,GLB_ALL);
  glbSwitchSystematics(GLB_ALL,GLB_ALL,GLB_ON);   
  fprintf(outfile1,"chi2 with statistics only: %g\n\n",chi2);

  glbDefineParams(input_errors,theta12*0.1,0,0,0,dm21*0.1,dm32/3);
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbDefineParams(central_values,theta12,theta13,theta23,deltacp,dm21,-dm32);
  glbSetDensityParams(central_values,1.0,GLB_ALL);
  glbSetCentralValues(central_values);
  glbSetInputErrors(input_errors);
  chi2=glbChiAll(central_values,minimum,GLB_ALL); 
  fprintf(outfile1,"chi2 at minimum: %g \n",chi2);
  fprintf(outfile1,"Position of minimum: theta12,theta13,theta23,delta,dm21,dm32,rho\n");
  glbPrintParams(outfile1,minimum);  

  /*---------------------------- Test 2 --------------------------*/
  
  printf("Test 2\n");
  double thetheta13,x,res1,res2;    
  for(x=-4;x<-2.0+0.001;x=x+2.0/50)
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
      
      fprintf(outfile2, "%g %g %g \n", x, res1, res2);
    }

  /*---------------------------- Test 3 --------------------------*/

  printf("Test 3\n");

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

        fprintf(outfile3, "%g %g %g \n",x,y,res);
    }
  }

  glbFreeParams(true_values);
  glbFreeParams(fit_values); 
  glbFreeParams(central_values);
  glbFreeParams(input_errors);
  glbFreeParams(minimum);
  glbFreeProjection(myprojection);
  glbFreeProjection(theta13_projection);
  /*---------------------------- Test 4 --------------------------*/

  printf("Test 4\n");

  /* Compute 1st edges of bars: systematics off, fit value of deltacp=0 */
  for(x=-7.0;x<-2.0+0.01;x=x+5.0/50)
  {
     res1=CalcNoSystematics(2.0e-3,x);
     res2=CalcNoSystematics(3.0e-3,x);
     fprintf(outfile4, "%g %g %g \n", x, res1, res2);
  }
  
  /* Compute 2nd edges of bars: systematics on, fit value of deltacp=0 */
  for(x=-7.0;x<-2.0+0.01;x+=5.0/50)
  {
     res1=CalcSystematics(2.0e-3,x);
     res2=CalcSystematics(3.0e-3,x);
     fprintf(outfile5, "%g %g %g \n", x, res1, res2);
  }
 
  /* Compute 3rd edges of bars: systematics on, projection onto s22th13 axis */
  glbDefineParams(original1,theta12,theta13,theta23,deltacp,dm21,2.0e-3);
  glbSetDensityParams(original1,1.0,GLB_ALL);
  glbDefineParams(original2,theta12,theta13,theta23,deltacp,dm21,3.0e-3);
  glbSetDensityParams(original2,1.0,GLB_ALL);
  for(x=-7.0;x<-2.0+0.01;x=x+5.0/50)
  {
     res1=CalcProjection(2.0e-3,x,original1);
     res2=CalcProjection(3.0e-3,x,original2);
     fprintf(outfile6, "%g %g %g \n", x, res1, res2);
  }

  /* Find sgn-degeneracies */
  FindDeg(degeneracy1,2.0e-3);
  FindDeg(degeneracy2,3.0e-3);
 
  /* Compute 4th edges of bars: systematics on, degeneracy, projection onto s22th13 axis */
  for(x=-7.0;x<-2.0+0.01;x=x+5.0/30)
  {
     res1=CalcProjection(2.0e-3,x,degeneracy1);
     res2=CalcProjection(3.0e-3,x,degeneracy2);
     fprintf(outfile7, "%g %g %g \n", x, res1, res2);
  }

  glbFreeParams(original1);
  glbFreeParams(original2);
  glbFreeParams(degeneracy1);
  glbFreeParams(degeneracy2);

  return 0;
}
