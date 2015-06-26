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

#include "myio.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static FILE* THEFILE;


void mioInitOutput(char* filename)
{
 if(strlen(filename)==0) THEFILE=stdout;
 else 
 {
   THEFILE=fopen(filename, "w");
   if (!THEFILE)
   {
     printf("File cannot be opened!\n");
     THEFILE=stdout;
   }
 }
 if(THEFILE!=stdout) fflush(THEFILE);
}

void mioCloseOutput()
{
 if(THEFILE!=stdout) fclose(THEFILE);
}

void mioAddToOutput(double n1,double n2)
{
 fprintf(THEFILE,"%g %g\n",n1,n2);
 if(THEFILE!=stdout) fflush(THEFILE);
}


