/* 
  File: main.c
  Author: Anthony R. Cassandra
  August, 1998

  ***
  Copyright 1994-1997, Brown University
  Copyright 1998, Anthony R. Cassandra

                           All Rights Reserved
                           
  Permission to use, copy, modify, and distribute this software and its
  documentation for any purpose other than its incorporation into a
  commercial product is hereby granted without fee, provided that the
  above copyright notice appear in all copies and that both that
  copyright notice and this permission notice appear in supporting
  documentation.
  
  ANTHONY CASSANDRA DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
  INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
  PARTICULAR PURPOSE.  IN NO EVENT SHALL ANTHONY CASSANDRA BE LIABLE FOR
  ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
  *****

  This file contains the main routine for the pomdp-solve program.
  
  The command line arguments are shown by running:

        pomdp-solve -h
*/
#define MAIN_C

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <mdp.h>

#include "global.h"
#include "signal-handler.h"
#include "cmd-line.h"
#include "pomdp.h"
#include "alpha.h"
#include "stats.h"
#include "random.h"
#include "policy-graph.h"
#include "projection.h"
#include "enumeration.h"
#include "linear-support.h"
#include "two-pass.h"
#include "witness.h"
#include "inc-prune.h"
#include "pomdp-solve.h"

/**********************************************************************/
void showUsagePomdpSolve( FILE *file ) {
/*
  Displays usage information in response to the '-h' or -'help'
  command line options. 
*/

  fprintf( file, "Usage: %s\n", gExecutableName );

  showUsageParams( file );

  /*******************************/
  /* Verbose flags               */
  /*******************************/
  fprintf( file, "Debug options:\n" );

  showUsageEnumType( file,
                     CMD_ARG_VERBOSE,
                     NUM_VERBOSE_MODES,
                     verbose_mode_str );

}  /* showUsagePomdpSolve */
/**********************************************************************/
PomdpSolveParams parseCmdLine( int argc, char **argv ) {
/*
  Parses the pomdp-solve command line, setting variablkes for the
  different settings.  The actual actions taken because of these
  settings are done later.

  Returns '1' if successful and '0' if there was a problem. 
*/
  PomdpSolveParams param;
  char str[MAX_MSG_LENGTH];
  int *mark_args;
  int i;
  int cmd_line_valid = TRUE;
  
   /****************/
  /* Copy the name of the executable into the global string.  This
     will be used for usage information so we need to do it right
     away in case something goes wrong reading the command line args. */ 
  strcpy( gExecutableName, argv[0] );
  
   /****************/
  /* We want to allow the different algorithms to parse their own
     arguments, but globally we want to make sure that all the
     command line arguments are valid.  To do this, we create a
     vector with an element for each command line argument.  We
     initialize it to zero and send it to all the individual command
     line parsing routines.  These routines will increment the
     counter for each option they actually use.  When we are all done,
     we can check to make sure everthing has been used and flag those
     unreckognized options.  Nothing is worse than putting a typo on
     the command line, thinking you have selected the right option,
     but in reality the program is doing something different. */
  mark_args = (int *) calloc( argc, sizeof( int ));
  
   /****************/
  /* First let's look for the '-h' or '-help' options, and if they
     are there just say that we failed to parse the command line
     (i.e., return FALSE), since this will trigger the routine to
     show the usage information. */
  
  if ( getFlagParam( argc, argv,  CMD_ARG_HELP_SHORT, mark_args )
       || getFlagParam( argc, argv, CMD_ARG_HELP_LONG, mark_args )) {
    showUsagePomdpSolve( stdout );
    exit ( -1 );
  }
  
   /****************/
   /* First see if a random number seed is given, and if it is set the
      seed to this value. */
   if( getStringParam( argc, argv, 
                       CMD_ARG_RAND_SEED,
                       mark_args, str )) 
     setRandomSeedFromString( str );

   /* Otherwise initialize the random number generator with
      psuedo-random seed. */
   else
     randomize();
   
   /****************/
   /* Set what verbose modes are specified.  The verbose argument is a
      comma separated list with *no* white-space. */

   /* Start by assumning all verbose modes are off. */
   for ( i = 0; i < NUM_VERBOSE_MODES; i++ )
     gVerbose[i] = FALSE;

   if( getStringParam( argc, argv, 
                       CMD_ARG_VERBOSE,
                       mark_args, str ))
     parseVerboseModes( str );

   /* Now we get all the POMDP solving specific parameters. */
   param = parseCmdLineParams( argc, argv, mark_args );

   /* Now we look for arguments that weren't parsed or which
      were parsed more than once. The first argument is just
      the program name, so start the loop at 1. */

   for ( i = 1; i < argc; i++ )
     if ( mark_args[i] != 1 ) {
       
       fprintf( param->report_file, 
                "** Error: cmd line option problem for '%s'\n",
                argv[i] );
       cmd_line_valid = FALSE;
     }
   
   free( mark_args );

   /* The one required command line option is the POMDP filename.
      Without it, it makes no sense to continue. */
   if ( param->param_filename[0] == NULL_CHAR ) {
     fprintf( param->report_file, 
              "** Error: no POMDP file specified ('-p' option)\n" );
     cmd_line_valid = FALSE;
   }

   /* If we had a problem parsing the command line, then show the
      help message and bag out. */
   if ( ! cmd_line_valid ) {

     fprintf
       ( param->report_file, 
         "Use %s or %s for the full set of command line options.\n",
         CMD_ARG_HELP_SHORT, CMD_ARG_HELP_LONG );
     exit( -1 );

   } /* if bad command line */  

   return ( param );

}  /* parseCmdLine */
/**********************************************************************/
int main( int argc, char **argv ) {
  PomdpSolveParams param;
  
  /* First thing is to parse command line to check for validity and to
     extract all the things we need. The program will exit and/or
     print a message indicating problems with the command line
     parsing, so we need not take any action here. */ 
  param = parseCmdLine( argc, argv );

  /* Do any file reading and memory allocations that are required. */
  initPomdpSolve( param );
  
  /* To document what is actually being executed we dump out
     all the parameters of the execution. */
  showPomdpSolveParams( param );

  /* And away we go... */
  solvePomdp( param );
  
  /* Let's make all nice afterwards for all other things. Note that
     this deallocates the space for 'param'. */
  cleanUpPomdpSolve( param );

  return( 1 );

} /* main */
/**********************************************************************/

