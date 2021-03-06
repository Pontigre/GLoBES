/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */
#include "myio.h"	     /* Simple output tools */

int main(int argc, char *argv[])
{
  /* Initialize libglobes */
  glbInit(argv[0]);

  glbInitExperiment("NFstandard.glb",&glb_experiment_list[0],&glb_num_of_exps);

  glbSelectMinimizer(GLB_MIN_POWELL); // Use experimental minimizer to speed up things

  double theta12 = asin(sqrt(0.3));
  double theta13 = asin(sqrt(0.001))/2.0;
  double theta23 = M_PI/4.0;
  double deltacp = M_PI/4.0;
  double sdm = 7.0e-5;
  double ldm = 2.0e-3;
  
  /* Parameter is output file name; if empty string, write to screen: */
  mioInitOutput("tut4.dat"); 

  glb_params fit_values = glbAllocParams();
  glb_params true_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_params foundmin = glbAllocParams();
  
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(true_values,1.0,GLB_ALL);  		 
  glbDefineParams(input_errors,theta12*0.1,0,0,0,sdm*0.1,0);     /* 10% external error for solar parameters */ 
  glbSetDensityParams(input_errors,0.05,GLB_ALL);  		 /* 5% matter density uncertainty */

  glbSetOscillationParameters(true_values);
  glbSetRates();
 
  glbSetCentralValues(true_values);
  glbSetInputErrors(input_errors);
 
  double res;

  /* scan parameter space as function of true x = log sin^2 2 theta_13 */
  double x;  
  double minx = -3.5;
  double maxx = -2.0;
  double steps = 18.0;

  glbCopyParams(true_values,fit_values);

  for(x=minx;x<maxx+0.0001;x+=(maxx-minx)/steps)
  {
     glbSetOscParams(fit_values,asin(sqrt(pow(10,x))/2),GLB_THETA_13);

     res=glbChiTheta13(fit_values,foundmin,GLB_ALL);

     printf("\nx=%g, Minimum found at %g:\n",x,res);
     glbPrintParams(stdout,foundmin);
 
     /* Write to file */
     mioAddToOutput(x,res);
  }
 
  /* Close output stream: */
  mioCloseOutput();
    
  glbFreeParams(fit_values);
  glbFreeParams(true_values);
  glbFreeParams(input_errors);
  glbFreeParams(foundmin);
  
 
  exit(0);
}

