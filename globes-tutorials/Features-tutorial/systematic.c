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

/* safe system calls */

FILE *my_fopen(const char *path, const char *mode)
{
  FILE *fp;
  fp=fopen(path,mode);
  if(fp==NULL) {fprintf(stderr,"FATAL: Could not open file '%s'!\n",path);exit(1);}
  return fp;
}

#define fopen my_fopen

/* "True" oscillation parameters */
double theta12,theta13,theta23,deltacp,sdm,ldm;

/* GLoBES parameter structures */
glb_params true_values;
glb_params test_values;
glb_params input_errors;

#define MAX_SYS 200
double sys_errors[MAX_SYS];       /* Uncertainties of systematics params */
double sys_startval[MAX_SYS];     /* Starting values for systematics minimizer */
double sigma_binbin = 0.0;        /* Bin-to-bin error */
double step_size=0.01;            /* step size in log sin^22theta_{13} */


#define EXP_FAR  0
#define EXP_NEAR 1


/***************************************************************************
 *                        H E L P E R   F U N C T I O N S                  *
 ***************************************************************************/

/* Minimum of two numbers */
inline double min(double x, double y)
{
  if (x < y)
    return x;
  else
    return y;
}

/* Square of real number */
inline double square(double x)
{
  return x*x;
}

/* Gauss likelihood (this is sufficient here due to the large event numbers
 * in a reactor experiment; for other setups, one should use Poisson statistics) */
inline double likelihood(double true_rate, double fit_rate, double sqr_sigma)
{
  if (sqr_sigma > 0)
    return square(true_rate - fit_rate) / sqr_sigma;
  else
    return 0.0;
}



/***************************************************************************
 *                      C H I ^ 2   F U N C T I O N S                      *
 ***************************************************************************/

/***************************************************************************
 * Calculate chi^2 for Double Chooz, including the following systematical  *
 * errors:                                                                 *
 *   x[0]: Flux normalization of reactor                                   *
 *   x[1]: Fiducial mass error - far detector                              *
 *   x[2]: Fiducial mass error - near detector                             *
 *   x[3]: Energy calibration error - far detector                         *
 *   x[4]: Energy calibration error - near detector                        *
 ***************************************************************************/
double chiDCNorm(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data)
{
  int n_bins = glbGetNumberOfBins(EXP_FAR);
  double *true_rates_N = glbGetRuleRatePtr(EXP_NEAR, 0);
  double *true_rates_F = glbGetRuleRatePtr(EXP_FAR, 0);
  double signal_fit_rates_N[n_bins];
  double signal_fit_rates_F[n_bins];
  double signal_norm_N, signal_norm_F;
  int ew_low, ew_high;
  double emin, emax;
  double fit_rate;
  double chi2 = 0.0;
  int i;

  /* Request simulated energy interval and analysis energy window */
  glbGetEminEmax(exp, &emin, &emax);
  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);

  /* Apply energy calibration error */
  glbShiftEnergyScale(x[3], glbGetSignalFitRatePtr(EXP_FAR, 0),
                      signal_fit_rates_F, n_bins, emin, emax);
  glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_NEAR, 0),
                      signal_fit_rates_N, n_bins, emin, emax);

  /* Loop over all bins in energy window */
  signal_norm_F = 1.0 + x[0] + x[1];
  signal_norm_N = 1.0 + x[0] + x[2];
  for (i=ew_low; i <= ew_high; i++)
  {
    /* Statistical part of chi^2 for far detector
     * Normalization is affected by flux error x[0] and fiducial mass error x[1] */
    fit_rate  = signal_norm_F * signal_fit_rates_F[i];
    chi2 += likelihood(true_rates_F[i], fit_rate, true_rates_F[i]);

    /* Statistical part of chi^2 for near detector
     * Normalization is affected by flux error x[0] and fiducial mass error x[2] */
    fit_rate  = signal_norm_N * signal_fit_rates_N[i];
    chi2 += likelihood(true_rates_N[i], fit_rate, true_rates_N[i]);
  }

  /* Systematical part of chi^2 (= priors) */
  for (i=0; i < n_params; i++)
    chi2 += square(x[i] / errors[i]);

  /* Save the systematics parameters as starting values for the next step */
  for (i=0; i < n_params; i++)
    sys_startval[i] = x[i];

  return chi2;
}  


/***************************************************************************
 * Calculate chi^2 for Double Chooz, including the following systematical  *
 * errors:                                                                 *
 *   x[0]: Flux normalization of reactor                                   *
 *   x[1]: Fiducial mass error - far detector                              *
 *   x[2]: Fiducial mass error - near detector                             *
 *   x[3]: Energy calibration error - far detector                         *
 *   x[4]: Energy calibration error - near detector                        *
 *   x[5]...x[5+nBins-1]: Spectral error                                   *
 * and the bin-to-bin error sigma_binbin. The calculation is based on      *
 * Eq. (9) of hep-ph/0303232.                                              *
 ***************************************************************************/
double chiDCSpectral(int exp, int rule, int n_params, double *x, double *errors,
                     void *user_data)
{
  int n_bins = glbGetNumberOfBins(EXP_FAR);
  double *true_rates_N = glbGetRuleRatePtr(EXP_NEAR, 0);
  double *true_rates_F = glbGetRuleRatePtr(EXP_FAR, 0);
  double signal_fit_rates_N[n_bins];
  double signal_fit_rates_F[n_bins];
  double signal_norm_N, signal_norm_F;
  int ew_low, ew_high;
  double emin, emax;
  double fit_rate;
  double chi2 = 0.0;
  int i;

  /* Request simulated energy interval and analysis energy window */
  glbGetEminEmax(exp, &emin, &emax);
  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);

  /* Apply energy calibration error */
  glbShiftEnergyScale(x[3], glbGetSignalFitRatePtr(EXP_FAR, 0),
                      signal_fit_rates_F, n_bins, emin, emax);
  glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_NEAR, 0),
                      signal_fit_rates_N, n_bins, emin, emax);

  /* Loop over all bins in energy window */
  signal_norm_F = 1.0 + x[0] + x[1];
  signal_norm_N = 1.0 + x[0] + x[2];
  for (i=ew_low; i <= ew_high; i++)
  {
    /* Statistical part of chi^2 for far detector
     * Normalization is affected by flux error x[0], fiducial mass error x[1],
     * spectral perturbations x[5+i], and bin-to-bin error sigma_binbin */
    fit_rate  = (signal_norm_F + x[5+i]) * signal_fit_rates_F[i];
    chi2 += likelihood( true_rates_F[i], fit_rate, 
                    true_rates_F[i] * (1.0 + true_rates_F[i]*square(sigma_binbin)) );

    /* Statistical part of chi^2 for near detector
     * Normalization is affected by flux error x[0], fiducial mass error x[2],
     * spectral perturbations x[5+i], and bin-to-bin error sigma_binbin */
    fit_rate  = (signal_norm_N + x[5+i]) * signal_fit_rates_N[i];
    chi2 += likelihood( true_rates_N[i], fit_rate, 
                    true_rates_N[i] * (1.0 + true_rates_N[i]*square(sigma_binbin)) );
  }

  /* Systematical part of chi^2 (= priors) */
  for (i=0; i < n_params; i++)
    chi2 += square(x[i] / errors[i]);

  /* Save the systematics parameters as starting values for the next step */
  for (i=0; i < n_params; i++)
    sys_startval[i] = x[i];

  return chi2;
}  



/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/

int main(int argc, char *argv[])
{ 
  double *old_sys_errors = NULL;      /* Temp. pointer to systematical error array */
  int sys_dim;                        /* Abbrv. for number of systematical errors */
  int n_bins=62;                      /* Number of bins */
  int i;
  
  /* Initialization */
  for (i=0; i < MAX_SYS; i++)
    sys_startval[i] = 0.0;

  /* Set standard oscillation parameters (cf. hep-ph/0405172v5) */
  theta12 = asin(sqrt(0.3));
  theta13 = 0.0;
  theta23 = M_PI/4;
  deltacp = M_PI/2;
  sdm = 7.9e-5;
  ldm = 2.6e-3;

  glbInit(argv[0]);                    /* Initialize GLoBES and define chi^2 functions */
  glbDefineChiFunction(&chiDCNorm,     5,        "chiDCNorm",     NULL);
  glbDefineChiFunction(&chiDCSpectral, 5+n_bins, "chiDCSpectral", NULL);
 
  /* Load 2 experiments: DC far (#0) and near (#1) detectors */
  glbClearExperimentList();
  glbInitExperiment("D-Chooz_far.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment("D-Chooz_near.glb", &glb_experiment_list[0], &glb_num_of_exps);
  if (glbGetNumberOfBins(EXP_FAR) != n_bins || glbGetNumberOfBins(EXP_NEAR) != n_bins)
  {
    printf("ERROR: Number of bins changed in AEDL file, but not in C code (or vice-versa).\n");
    return -1;
  }
  else
    n_bins = glbGetNumberOfBins(EXP_FAR);


  /* Initialize parameter vectors */
  true_values  = glbAllocParams();
  test_values  = glbAllocParams();
  input_errors = glbAllocParams();
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  
  glbSetDensityParams(test_values,1.0,GLB_ALL);
  glbDefineParams(input_errors, 0.1*theta12, 0, 0.15*theta23, 0, 0.05*sdm, 0.05*ldm);
  glbSetDensityParams(input_errors, 0.05, GLB_ALL);
  glbSetOscillationParameters(true_values);
  glbSetInputErrors(input_errors);

  /* Compute chi^2 as a function of sin^22theta */


  /* Calculate "true" event rates */
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  
    glbSetDensityParams(test_values,1.0,GLB_ALL);
    glbDefineParams(input_errors, 0.1*theta12, 0, 0.15*theta23, 0, 0.05*sdm, 0.05*ldm);
    glbSetDensityParams(input_errors, 0.05, GLB_ALL);
    glbSetInputErrors(input_errors);   

    glbSetOscillationParameters(true_values);
    glbSetRates();

    FILE *fp=fopen("sys-data0","w");

    double thetheta13, chi2, x ;
    /* First without sytematics */
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF);
    fprintf(fp,"# no systematics\n"); 
    for(x=-3;x<=-1;x+=step_size)
      {
	
	/* Set vector of test values */
	thetheta13 = asin(sqrt(pow(10,x)))/2;
	glbSetOscParams(test_values, thetheta13, GLB_THETA_13);
	
	/* Set starting values for systematics minimiyer to the coordinates of
	 * minimum in the last iteration. This accelerates the minimization and
	 * prevents convergence problems. */
	glbSetSysStartingValuesList(EXP_FAR, 0, GLB_ON, sys_startval);
	
	/* Compute Chi^2 for all loaded experiments and all rules
	 * Correlations are unimportant in reactor experiments, so glbChiSys is sufficient */
	chi2 = glbChiSys(test_values, GLB_ALL, GLB_ALL);
	fprintf(fp,"%f\t%f\n",x,chi2);
      }

    fclose(fp);

    /* Calculate sensitivity curve with normalization and energy calibration errors,
     * as defined in the AEDL files */

    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON);
    glbSetOscillationParameters(true_values);
    glbSetRates();


    fp=fopen("sys-data1","w");
    fprintf(fp,"# normalization & energy calibration as in the AEDL files\n"); 
    for(x=-3;x<=-1;x+=step_size)
      {
	
	/* Set vector of test values */
	thetheta13 = asin(sqrt(pow(10,x)))/2;
	glbSetOscParams(test_values, thetheta13, GLB_THETA_13);
	
	/* Set starting values for systematics minimizer to the coordinates of
	 * minimum in the last iteration. This accelerates the minimization and
	 * prevents convergence problems. */
	glbSetSysStartingValuesList(EXP_FAR, 0, GLB_ON, sys_startval);
	
	/* Compute Chi^2 for all loaded experiments and all rules
	 * Correlations are unimportant in reactor experiments, so glbChiSys is sufficient */
	chi2 = glbChiSys(test_values, GLB_ALL, GLB_ALL);
	fprintf(fp,"%f\t%f\n",x,chi2);
      }


    fclose(fp);
    fp=fopen("sys-data2","w");

    /* Calculate sensitivity curve with the above + spectral error
     * Since chiDCSpectral computes the complete chi^2 for the whole problem, it
     * must be called only for ONE of the two experiments */

    old_sys_errors = glbGetSysErrorsListPtr(EXP_FAR, 0, GLB_ON);   /* Fill error array */
    sys_dim        = glbGetSysDimInExperiment(EXP_FAR, 0, GLB_ON);
    for (i=0; i < sys_dim; i++)         /* Normalization and energy calibration errors */
      sys_errors[i] = old_sys_errors[i];
     for (i=sys_dim; i < sys_dim + n_bins; i++)
       sys_errors[i] = 0.02;                                          /* Spectral error */
    sigma_binbin = 0.0;                          /* No bin-to-bin error for the moment */
    glbSetChiFunction(EXP_FAR, 0, GLB_ON, "chiDCSpectral", sys_errors);
    glbSetChiFunction(EXP_NEAR, 0, GLB_ON, "chiZero", sys_errors);
    
    glbSetOscillationParameters(true_values);
    glbSetRates();

    fprintf(fp,"#as above + spectral error\n"); 
    for(x=-3;x<=-1;x+=step_size)
      {
	
	/* Set vector of test values */
	thetheta13 = asin(sqrt(pow(10,x)))/2;
	glbSetOscParams(test_values, thetheta13, GLB_THETA_13);
	
	/* Set starting values for systematics minimiyer to the coordinates of
	 * minimum in the last iteration. This accelerates the minimization and
	 * prevents convergence problems. */
	glbSetSysStartingValuesList(EXP_FAR, 0, GLB_ON, sys_startval);
	
	/* Compute Chi^2 for all loaded experiments and all rules
	 * Correlations are unimportant in reactor experiments, so glbChiSys is sufficient */
	chi2 = glbChiSys(test_values, GLB_ALL, GLB_ALL);
	fprintf(fp,"%f\t%f\n",x,chi2);
      }

    fclose(fp);
    fp=fopen("sys-data3","w");


    old_sys_errors = glbGetSysErrorsListPtr(EXP_FAR, 0, GLB_ON);   /* Fill error array */
    sys_dim        = glbGetSysDimInExperiment(EXP_FAR, 0, GLB_ON);
    for (i=0; i < sys_dim; i++)         /* Normalization and energy calibration errors */
      sys_errors[i] = old_sys_errors[i];
    for (i=sys_dim; i < sys_dim + n_bins; i++)
      sys_errors[i] = 0.02;                                          /* Spectral error */
    sigma_binbin = 0.02;                          /* No bin-to-bin error for the moment */
    glbSetChiFunction(EXP_FAR, 0, GLB_ON, "chiDCSpectral", sys_errors);
    glbSetChiFunction(EXP_NEAR, 0, GLB_ON, "chiZero", sys_errors);
    
    glbSetOscillationParameters(true_values);
    glbSetRates();


    fprintf(fp,"# as above + bin-to-bin error\n"); 
    for(x=-3;x<=-1;x+=step_size)
      {
	
	/* Set vector of test values */
	thetheta13 = asin(sqrt(pow(10,x)))/2;
	glbSetOscParams(test_values, thetheta13, GLB_THETA_13);
	
	/* Set starting values for systematics minimiyer to the coordinates of
	 * minimum in the last iteration. This accelerates the minimization and
	 * prevents convergence problems. */
	glbSetSysStartingValuesList(EXP_FAR, 0, GLB_ON, sys_startval);
	
	/* Compute Chi^2 for all loaded experiments and all rules
	 * Correlations are unimportant in reactor experiments, so glbChiSys is sufficient */
	chi2 = glbChiSys(test_values, GLB_ALL, GLB_ALL);
	fprintf(fp,"%f\t%f\n",x,chi2);
      }

    fclose(fp);
  
  /* Clean up */
  glbFreeParams(true_values);
  glbFreeParams(test_values); 
  glbFreeParams(input_errors);
  
  return 0;  
}


