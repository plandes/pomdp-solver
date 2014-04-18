######################################################################
#
#  File: Makefile
#  Author: A. R. Cassandra
#  July, 1998
#
#  *****
#  Copyright 1994-1997, Brown University
#  Copyright 1998, Anthony R. Cassandra
#
#                           All Rights Reserved
#                           
#  Permission to use, copy, modify, and distribute this software and its
#  documentation for any purpose other than its incorporation into a
#  commercial product is hereby granted without fee, provided that the
#  above copyright notice appear in all copies and that both that
#  copyright notice and this permission notice appear in supporting
#  documentation.
#  
#  ANTHONY CASSANDRA DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
#  INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
#  PARTICULAR PURPOSE.  IN NO EVENT SHALL ANTHONY CASSANDRA BE LIABLE FOR
#  ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
#  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
#  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
#  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#  *****
#
######################################################################

# See the comments that start with 4 pound (#) signs for things that
# may need to change for system-specific settings.

######################################################################
# Executable name and location settings

# NOTE: If you do not use the GNU compiler 'gcc' then define this
# constant to be -DNON_GNU, since it will alleviate some portability
# problems.  If you use the GNU compiler, then leave this undefined or
# defined to be empty.
#
#GNU_DEFINES = -DNON_GNU
GNU_DEFINES = 

# The compiler and some flags for dependencies
ifneq ($(GNU_DEFINES),-DNON_GNU)
CC = gcc
DEPEND_FLAGS = -E -M
else
CC = cc
DEPEND_FLAGS = -xM
endif

# For scanning
LEX = lex

# For parsing
YACC = yacc -v -d

# Removing files
RM = /bin/rm -f

# Copying files
CP = /bin/cp

# Archiver name and flags
AR = ar -crv

######################################################################
# OS specific flags

# Flags specific to the machine architecture.  The CPLEX tech-support
# people tell me that I should use SYSSOLARIS on the Suns.  I am not
# sure it matters.

###################################
# OS Setting

# Make a variable depending on the operating system type.

SYS_FLAGS = -DSYSUNKNOWN

ifeq ($(OSTYPE),solaris)
SYS_FLAGS = -DSYSSOLARIS
endif

ifeq ($(OSTYPE),linux)
SYS_FLAGS = -DSYSLINUX
endif

######################################################################
# LP solver settings

# Assume that the presence of the environment variable CPLEXLICENSE
# means CPLEX is present.

ifndef CPLEXLICENSE

# If no CPLEX then use lp_solve
LP_DEFINES = 
CPLEX_INCLUDE_PATH =
CPLEX_LIB_LD_PATH = 
CPLEX_LIBS = 
CPLEX_EXTRA_OBJS = 

endif

# You probably need to adjust the CPLEX base location for your system.
# You can set it as an environment variable or set it here.

ifndef CPLEX_HOME
CPLEX_HOME = /pro/cplex
endif

ifeq ($(CPLEX_VERSION),3)

# For CPLEX 3.0
LP_DEFINES = -DHAVE_CPLEX 
CPLEX_INCLUDE_PATH = -I$(CPLEX_HOME)/include
CPLEX_LIB_LD_PATH =
CPLEX_LIBS = $(CPLEX_HOME)/lib/cplex.a
CPLEX_EXTRA_OBJS = 

endif

ifeq ($(CPLEX_VERSION),4)

# For CPLEX 4.0
LP_DEFINES = -DHAVE_CPLEX
CPLEX_INCLUDE_PATH = -I$(CPLEX_HOME)
CPLEX_LIB_LD_PATH = -L$(CPLEX_HOME)
CPLEX_LIBS = -lcplex
CPLEX_EXTRA_OBJS =  -DCPX_PROTOTYPE_ANSI $(CPLEX_HOME)/oldcplex.c

endif


###################################
# LP_SOLVE

# We compile lp_solve even if we have CPLEX, so here are the constants
# we need.   Note that lp_solve will be compiled but not used.

LP_SOLVE_LIB_PATH = ./lp_solve
LP_SOLVE_LIB_NAME = lpk
LP_SOLVE_INCLUDE_PATH = -I$(LP_SOLVE_LIB_PATH)
LP_SOLVE_LIB_LD_PATH = -L$(LP_SOLVE_LIB_PATH)
LP_SOLVE_LIBS = -llpk

######################################################################
# Program options
#
# Add these to the POMDP_SOLVE_DEFINES if wanted.
#
#   -DUSE_DENSE_LPS  
#       Force the program to represent the LPs densely.
#       (see region.h)
#
#   -DUSE_OLD_LP_FORMULATION 
#       Force the use of the old-style region LP constraints.
#       (see region.h)
# 
#   -DDISABLE_MEMORY_LIMIT
#   -DDISABLE_TIME_LIMIT
#   -DDISABLE_SIGNAL_HANDLING

MDP_LIB_PATH = ./mdp
MDP_LIB_NAME = mdp
MDP_INCLUDES = -I$(MDP_LIB_PATH)
MDP_LIB_LD_PATH = -L$(MDP_LIB_PATH)
MDP_LIBS = -l$(MDP_LIB_NAME)

######################################################################
# Compiler - compile time flags

#PROFILE_FLAGS = -pg
#OPTIMIZE_FLAGS = -O3

CFLAGS = \
	$(PROFILE_FLAGS) \
	$(OPTIMIZE_FLAGS) \
	$(SYS_FLAGS) \
	$(LP_DEFINES) \
	$(GNU_DEFINES)

# Extra load and include paths
INCLUDES = $(MDP_INCLUDES) $(LP_SOLVE_INCLUDE_PATH) $(CPLEX_INCLUDE_PATH)

###################################
# Compiler load-time flags

SYS_LIBS = -lm

LIB_PATH = $(MDP_LIB_LD_PATH) $(LP_SOLVE_LIB_LD_PATH) $(CPLEX_LIB_LD_PATH)

CLDLIBS  = $(MDP_LIBS) $(LP_SOLVE_LIBS) $(CPLEX_LIBS) $(SYS_LIBS)

############################################################
# Special debugging libraries 
#

# For the electric fence libraries.
#
#DEBUGLIBS=./tools/ElectricFence-2.0.5/libefence.a

DEBUG_FLAGS = -g -pedantic

# For the dmalloc libraries
#
#DEBUGLIBS= -ldmalloc
#DEBUG_FLAGS = -g -DUSE_DMALLOC

############################################################
# Main program binary and object files
POMDP_SOLVE_BIN = pomdp-solve
POMDP_SOLVE_OBJS = \
	global.o \
	timing.o \
	random.o \
	pomdp.o \
	alpha.o \
	stats.o \
	params.o \
	projection.o \
	lp-interface.o \
	region.o \
	parsimonious.o \
	neighbor.o \
	common.o \
	witness.o \
	cross-sum.o \
	enumeration.o \
	inc-prune.o \
	two-pass.o \
	vertex-enum.o \
	linear-support.o \
	policy-graph.o \
	cmd-line.o \
	signal-handler.o \
	pomdp-solve.o \
	main.o

############################################################
# Test program binary and object files
AUTO_TEST_BIN = auto-test
AUTO_TEST_OBJS = \
	timing.o \
	random.o \
	global.o \
	pomdp.o \
	alpha.o \
	stats.o \
	cmd-line.o \
	projection.o \
	lp-interface.o \
	region.o \
	parsimonious.o \
	neighbor.o \
	common.o \
	enumeration.o \
	cross-sum.o \
	inc-prune.o \
	witness.o \
	two-pass.o \
	vertex-enum.o \
	linear-support.o \
	policy-graph.o \
	params.o \
	signal-handler.o \
	pomdp-solve.o \
	auto-test.o

############################################################

# Main executable
$(POMDP_SOLVE_BIN): libs $(POMDP_SOLVE_OBJS) 
	$(CC) -o $(POMDP_SOLVE_BIN) $(CLDFLAGS) $(LIB_PATH) $(POMDP_SOLVE_OBJS) $(CPLEX_EXTRA_OBJS) $(DEBUGLIBS) $(CLDLIBS)

# When we want to experiment with variations without changing the real
# executable (i.e., if some experiments are running.)
pomdp-solve-debug: $(POMDP_SOLVE_OBJS)
	$(CC) -o pomdp-solve-debug $(CLDFLAGS) $(LIB_PATH) $(POMDP_SOLVE_OBJS) $(CPLEX_EXTRA_OBJS) $(DEBUGLIBS) $(CLDLIBS)

# Libraries that are needed
$(MDP_LIB_PATH)/lib$(MDP_LIB_NAME).a:
	$(MAKE) -C $(MDP_LIB_PATH) lib$(MDP_LIB_NAME).a
$(LP_SOLVE_LIB_PATH)/lib$(LP_SOLVE_LIB_NAME).a:
	$(MAKE) -C $(LP_SOLVE_LIB_PATH) lib$(LP_SOLVE_LIB_NAME).a
libs:	$(MDP_LIB_PATH)/lib$(MDP_LIB_NAME).a \
		$(LP_SOLVE_LIB_PATH)/lib$(LP_SOLVE_LIB_NAME).a

# Automatic test program
$(AUTO_TEST_BIN):	libs $(AUTO_TEST_OBJS) 
	$(CC) -o $(AUTO_TEST_BIN) $(CLDFLAGS) $(LIB_PATH) $(AUTO_TEST_OBJS) $(CPLEX_EXTRA_OBJS) $(DEBUGLIBS) $(CLDLIBS)
	$(MAKE) -C lp_solve

test:	$(AUTO_TEST_BIN)
	$(AUTO_TEST_BIN)

# Everything
all:	$(POMDP_SOLVE_BIN) $(AUTO_TEST_BIN) test

######################################################################
# Do this to clean up before a fresh make

clean:
	$(MAKE) -C $(MDP_LIB_PATH) clean
	$(MAKE) -C $(LP_SOLVE_LIB_PATH) clean
	$(RM) *.o 

clean-local:
	$(RM) *.o 

######################################################################
# Specifies the rules

.SUFFIXES: .o .c .y .l .h

%.o:	%.c
	$(CC) $(POMDP_SOLVE_DEFINES) $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDES) -c $<

######################################################################
#  Specifies the dependencies

# Name of the dependency file
DEPEND = Makefile.depend

dep:
	$(MAKE) -C $(MDP_LIB_PATH) dep
	$(CC) $(DEPEND_FLAGS) $(INCLUDES) *.[c] > $(DEPEND)

# See if a Dependency file exists and if so use it
ifeq ($(DEPEND),$(wildcard $(DEPEND)))
include $(DEPEND)
endif

######################################################################


