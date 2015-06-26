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
 *   Hands-on session: Simulation of Accelerator neutrino experiments      *
 ***************************************************************************
 *   Part 1: Confidence regions in the th23-dm31 plane                     *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */

/* Output file */
char MYFILE[]="th23dm31.dat";
char AEDLFILE[]="T2K-tutorial.glb";
FILE *outfile = NULL;

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

  /* Define "true" oscillation parameters (cf. hep-ph/0405172v5) */
  double theta12 = asin(sqrt(0.3));
  double theta13 = asin(sqrt(0.0))/2.0;
  double theta23 = 45.0 * M_PI/180.0;
  double deltacp = 0.0;
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

  /* Compute simulated data */
  glbSetOscillationParameters(true_values);
  glbSetRates();

  /* Scan the th23-dm31 plane */
  double this_th23, this_ldm;
  double th23_lower = 35 * M_PI/180.0;
  double th23_upper = 56 * M_PI/180.0;
  double th23_steps = 20;
  double ldm_lower  = 2.4e-3;
  double ldm_upper  = 3.0e-3;
  double ldm_steps  = 20;
  double res;
  for(this_th23=th23_lower; this_th23 < th23_upper; this_th23 += (th23_upper-th23_lower)/th23_steps)
  {
    for(this_ldm=ldm_lower; this_ldm < ldm_upper; this_ldm += (ldm_upper-ldm_lower)/ldm_steps)
    {
      /* Set vector of test=fit values */
      glbSetOscParams(test_values, this_th23, GLB_THETA_23);
      glbSetOscParams(test_values, this_ldm, GLB_DM_31);

      /* Compute chi^2 including only systematical errors (no parameter correlations) */
      res = glbChiSys(test_values, GLB_ALL, GLB_ALL);
      fprintf(outfile, "%g %g %g\n", this_th23*180.0/M_PI, this_ldm, res);
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);
  
  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values);

  return 0;
}

