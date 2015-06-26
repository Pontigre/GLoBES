# GLoBES -- General LOng Baseline Experiment Simulator
# (C) 2002 - 2007,  The GLoBES Team
#
# GLoBES is mainly intended for academic purposes. Proper
# credit must be given if you use GLoBES or parts of it. Please
# read the section 'Credit' in the README file.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# The configure script should have replaced anything as required to
# obtain a working Makefile, which may serve as a template for own
# applications.
#
# This simple Makefile is for the GLoBES degeneracy tutorial
#
# Compile example files with ``make deg_tut_1'' etc.
#
# This Makefile assumes dynamic libraries, installed to either
# the default prefix /usr/local/ or to a user-defined directory
# called ${prefix}.
#
# For linking against a specific version of GLoBES, libglobes.so can be
# replaced by the respective library, such as libglobes.so.0.0.1
#
# For using static libraries, see the README file.

prefix = /home/winter/Projects/GLOBES/ModuleVersions
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

local_CFLAGS = -g -O4

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)

BIN = deg_tut_1 deg_tut_2 deg_tut_3 deg_tut_4 deg_tut_5 deg_tut_6 deg_tut_1_sol deg_tut_2_sol deg_tut_3_sol deg_tut_4_sol deg_tut_5_sol deg_tut_6_sol
OBJ = myio.o degtools.o deg_tut_1.o deg_tut_2.o deg_tut_3.o deg_tut_4.o deg_tut_5.o deg_tut_6.o deg_tut_1_sol.o deg_tut_2_sol.o deg_tut_3_sol.o deg_tut_4_sol.o deg_tut_5_sol.o deg_tut_6_sol.o

deg_tut_1: deg_tut_1.o myio.o 
	gcc deg_tut_1.o myio.o -o deg_tut_1 $(LDFLAGS) $(local_LDFLAGS)

deg_tut_1_sol: deg_tut_1_sol.o myio.o 
	gcc deg_tut_1_sol.o myio.o -o deg_tut_1_sol $(LDFLAGS) $(local_LDFLAGS)

deg_tut_2: deg_tut_2.o myio.o 
	gcc deg_tut_2.o myio.o -o deg_tut_2 $(LDFLAGS) $(local_LDFLAGS)

deg_tut_2_sol: deg_tut_2_sol.o myio.o 
	gcc deg_tut_2_sol.o myio.o -o deg_tut_2_sol $(LDFLAGS) $(local_LDFLAGS)

deg_tut_3: deg_tut_3.o myio.o 
	gcc deg_tut_3.o myio.o -o deg_tut_3 $(LDFLAGS) $(local_LDFLAGS)

deg_tut_3_sol: deg_tut_3_sol.o myio.o 
	gcc deg_tut_3_sol.o myio.o -o deg_tut_3_sol $(LDFLAGS) $(local_LDFLAGS)

deg_tut_4: deg_tut_4.o myio.o 
	gcc deg_tut_4.o myio.o -o deg_tut_4 $(LDFLAGS) $(local_LDFLAGS)

deg_tut_4_sol: deg_tut_4_sol.o myio.o 
	gcc deg_tut_4_sol.o myio.o -o deg_tut_4_sol $(LDFLAGS) $(local_LDFLAGS)

deg_tut_5: deg_tut_5.o myio.o 
	gcc deg_tut_5.o myio.o -o deg_tut_5 $(LDFLAGS) $(local_LDFLAGS)

deg_tut_5_sol: deg_tut_5_sol.o myio.o 
	gcc deg_tut_5_sol.o myio.o -o deg_tut_5_sol $(LDFLAGS) $(local_LDFLAGS)

deg_tut_6: deg_tut_6.o myio.o 
	gcc deg_tut_6.o myio.o -o deg_tut_6 $(LDFLAGS) $(local_LDFLAGS)

deg_tut_6_sol: deg_tut_6_sol.o myio.o 
	gcc deg_tut_6_sol.o myio.o -o deg_tut_6_sol $(LDFLAGS) $(local_LDFLAGS)

%.o : %.c
	gcc $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)
.PHONY: clean
clean:
	rm $(BIN) $(OBJ)
