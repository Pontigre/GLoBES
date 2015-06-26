/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2004,  The GLoBES Team
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

/***************************************************************************
 *   G L o B E S   -   P H Y S I C S   A N D   A P P L I C A T I O N S     *
 *                                                                         *
 *             24 - 26 January 2007, Heidelberg, Germany                   *
 *                                                                         *
 *              Hands-on session: Advanced Usage of AEDL                   *
 ***************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */

/* Output file */
char MYFILE[]="th13delta.dat";


/* ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' */
/*                                                                 */
/*            Please enter here the name of the glb-file           */
/*                                                                 */
/* ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' */
               
	       
  	             char AEDLFILE[]="NAME.glb";




FILE *outfile = NULL;

int main(int argc, char *argv[])
{ 
  /* Initialize libglobes */
  glbInit(argv[0]); 
  glbSelectMinimizer(GLB_MIN_POWELL);

  /* Initialize experiment */
  glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 
  glbInitExperiment("T2K_disappearance_only.glb",&glb_experiment_list[0],&glb_num_of_exps);


  /* Intitialize output */
  outfile = fopen(MYFILE, "w");
  if (outfile == NULL)
  {
    printf("Error opening output file.\n");
    return -1;
  }

  /* Define "true" oscillation parameters (cf. hep-ph/0405172v5) */
  double theta12 = asin(sqrt(0.3));
  double theta13 = asin(sqrt(0.1))/2.0;
  double theta23 = 45.0 * M_PI/180.0;
  double deltacp = 90.0 * M_PI/180.0;
  double sdm = 7.9e-5;
  double ldm = 2.6e-3;
  
  /* Define "true" oscillation parameter vector */
  glb_params true_values = glbAllocParams();
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
 
  /* Define initial guess for the fit values */ 
  glb_params test_values = glbAllocParams();
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  
  glbSetDensityParams(test_values,1.0,GLB_ALL);

  /* Define external input (1-sigma errors) on the parameters: 10% error
   * on the solar parameters, 5% on the matter density, all other parameters free.
   * External input is implemented as a prior of the form
   *    (fit_value - central_value)^2 / input_error^2
   */
  glb_params input_errors = glbAllocParams();
  glbDefineParams(input_errors, theta12*0.1, 0, 0, 0, sdm*0.1, 0);
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetInputErrors(input_errors);
  glbSetCentralValues(true_values);

  /* Define projection onto th13 and delta, marginalizing over
   * th23 and dm31. The solar parameters can be kept fixed to speed
   * up the calculation without introducing large errors. */
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
  double th13_lower  = 8*M_PI/180.;
  double th13_upper  = 11*M_PI/180.;
  double th13_steps  = 10;
  double delta_lower = 70*M_PI/180.;
  double delta_upper = 110*M_PI/180.;
  double delta_steps = 10;
  double res;
  for(this_th13=th13_lower; this_th13<=th13_upper; this_th13+=(th13_upper-th13_lower)/th13_steps)
  {
    for(this_delta=delta_lower; this_delta<=delta_upper; this_delta+=(delta_upper-delta_lower)/delta_steps)
    {
      /* Set vector of test=fit values */
      glbSetOscParams(test_values, this_th13, GLB_THETA_13);
      glbSetOscParams(test_values, this_delta, GLB_DELTA_CP);

      /* Compute chi^2 assuming the normal mass hierarchy in the fit */
      res = glbChiNP(test_values, NULL, GLB_ALL);
      fprintf(outfile, "%g %g %g\n", this_th13*180.0/M_PI, this_delta*180.0/M_PI, res);
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);
  
  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);
  glbFreeProjection(th13delta_projection);

  return 0;
}

