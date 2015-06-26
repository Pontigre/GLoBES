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

 /* SetExperiments(); */
  glbInitExperiment("T2HK.glb",&glb_experiment_list[0],&glb_num_of_exps); 

  double theta12 = asin(sqrt(0.3));
  double theta13 = asin(sqrt(0.1))/2.0;
  double theta23 = M_PI/4.0;
  double deltacp = 0.0;
  double sdm = 8.0e-5;
  double ldm = 2.5e-3;
  
  /* Parameter is output file name; if empty string, write to screen: */
  mioInitOutput("tut5.dat"); 

  glb_params central_values = glbAllocParams();
  glb_params fit_values = glbAllocParams();
  glb_params true_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_params foundmin = glbAllocParams();
  
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(true_values,1.0,GLB_ALL);  		 
  glbDefineParams(input_errors,theta12*0.1,0,0,0,sdm*0.1,0);     /* 10% external error for solar parameters */ 
  glbSetDensityParams(input_errors,0.05,GLB_ALL);  		 /* 5% matter density uncertainty */

  double x,res;
  for(x=0.45;x<0.55+0.0001;x+=0.01)
  {
   
    /* Compute chi2 at degeneracy for this choice of the simulated sin^2 theta_23: */
    glbSetOscParams(true_values,asin(sqrt(x)),GLB_THETA_23);
    
    /* Educated guess for fit values: */
    glbCopyParams(true_values,fit_values);
    glbSetOscParams(fit_values,M_PI/2.0-glbGetOscParams(true_values,GLB_THETA_23),GLB_THETA_23);
    glbCopyParams(fit_values,central_values);

    glbSetOscillationParameters(true_values);
    glbSetRates();
    glbSetInputErrors(input_errors);
    glbSetCentralValues(central_values);

    res=glbChiAll(fit_values,foundmin,GLB_ALL);

    printf("\nx=%g, Minimum found at %g:\n",x,res);
    glbPrintParams(stdout,foundmin);
    if(((glbGetOscParams(foundmin,GLB_THETA_23)<M_PI/4.0) && (glbGetOscParams(true_values,GLB_THETA_23)<M_PI/4.0)) || 
       ((glbGetOscParams(foundmin,GLB_THETA_23)>M_PI/4.0) && (glbGetOscParams(true_values,GLB_THETA_23)>M_PI/4.0)))
        printf("Wrong octant!\n");
 
    /* Write to file */
    mioAddToOutput(x,res);
  }

  /* Close output stream: */
  mioCloseOutput();
    
  glbFreeParams(fit_values);
  glbFreeParams(central_values);
  glbFreeParams(true_values);
  glbFreeParams(input_errors);
  glbFreeParams(foundmin);
 
  exit(0);
}

