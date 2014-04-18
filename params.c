/* 
  File: params.c 
  Author: Anthony R. Cassandra
  July, 1998

  *****
  Copyright 1994-1997, Brown University
  Copyright 1998, 1999 Anthony R. Cassandra

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
  
  Stuff to specify all POMDP solution parameters.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mdp.h>

#include "global.h"
#include "timing.h"
#include "random.h"
#include "pomdp.h"
#include "alpha.h"
#include "stats.h"
#include "cmd-line.h"
#include "lp-interface.h"
#include "enumeration.h"
#include "linear-support.h"
#include "two-pass.h"
#include "witness.h"
#include "inc-prune.h"
#include "params.h"

/* Strings for the various algorithms */
char *method_str[] = METHOD_STRINGS;

/* Strings for the various stopping criteria */
char *stop_criteria_str[] = STOP_CRITERIA_STRINGS;

/* Strings for the various stopping criteria */
char *purge_option_str[] = PURGE_OPTION_STRINGS;

/* Strings for the various incremental pruning variations. */
char *inc_prune_type_str[] = INC_PRUNE_TYPE_STRINGS;

/* Strings for the various value iterations variations. */
char *vi_variation_type_str[] = VI_VARIATION_TYPE_STRINGS;

/**********************************************************************/
PomdpSolveParams newPomdpSolveParams(  ) {
/*
  Creates the memory for the structure to hold the parameters used in
  solving a POMDP.  Also sets the fields to the default values.
*/
  PomdpSolveParams params;
  int i;

  params = malloc( sizeof( *params ));

  params->cur_epoch = 0;
  params->succinct = FALSE;
  params->stat_summary = FALSE;
  params->report_filename[0] = NULL_CHAR;
  params->report_file = stdout;
  params->param_filename[0] = NULL_CHAR;
  params->override_discount = -1.0;
  params->method = DEFAULT_METHOD;
  params->stop_criteria = DEFAULT_STOP_CRITERIA;
  params->horizon = DEFAULT_HORIZON;
  params->prefix_str[0] = NULL_CHAR;
  params->alpha_filename[0] = NULL_CHAR;
  params->pg_filename[0] = NULL_CHAR;
  params->initial_policy_filename[0] = NULL_CHAR;
  params->initial_policy = NULL;
  params->max_secs = 0;
  params->memory_limit = 0;
  params->save_all = FALSE;
  params->proj_purge = DEFAULT_PROJECTION_PURGE;
  params->use_witness_points = DEFAULT_USE_WITNESS_POINTS;
  params->backup_file[0] = NULL_CHAR;
  params->q_purge_option = DEFAULT_Q_PURGE_OPTION;
  params->domination_check = DEFAULT_DOMINATION_CHECK;
  params->alg_init_rand_points = DEFAULT_ALG_INIT_RAND_POINTS;
  params->prune_init_rand_points = DEFAULT_PRUNE_INIT_RAND_POINTS;

  params->prune_epsilon = DEFAULT_PRUNE_EPSILON; 
  params->epsilon = DEFAULT_EPSILON; 
  params->lp_epsilon = DEFAULT_LP_EPSILON; 

  params->weak_bound_delta = DEFAULT_WEAK_BOUND_DELTA;

  params->alpha_epsilon = DEFAULT_ALPHA_EPSILON; 
  params->vertex_epsilon = DEFAULT_VERTEX_EPSILON; 
  params->impossible_obs_epsilon = DEFAULT_IMPOSSIBLE_OBS_EPSILON;
  params->double_equality_precision = DEFAULT_DOUBLE_EQUALITY_PRECISION;
  
  /* Default value to use when considering whether to include a
     coefficient in a sparse representation. Note that we don't really
     want to tie this to the precision that is being used to solve the
     problem, because this value can change the problem being solved.
     Thus this should just be fixed for all time at the minimum
     precision. */
  params->sparse_epsilon = SMALLEST_PRECISION; 

  /* Place to hang statistics off of (optional) */
  params->stat = NULL;

  /****************************************/
  /****  Algorithm specific section  ******/
  /****************************************/

  params->ip_type = DEFAULT_INC_PRUNE_TYPE;
  params->enum_purge_option = DEFAULT_ENUM_PURGE_OPTION;

  /****************************************/
  /****  VI variation specific section  ***/
  /****************************************/

  params->vi_variation = DEFAULT_VI_VARIATION;
  params->starting_epsilon = DEFAULT_STARTING_EPSILON;
  params->ending_epsilon = DEFAULT_ENDING_EPSILON;
  params->epsilon_adjust_factor = DEFAULT_EPSILON_ADJUST_FACTOR;
  params->max_soln_size = DEFAULT_MAX_SOLN_SIZE;
  params->epoch_history_window_length = DEFAULT_HISTORY_WINDOW_LENGTH; 
  params->epoch_history_window_delta = DEFAULT_HISTORY_WINDO_DELTA;

  return ( params );

}  /* newPomdpSolveParams */
/**********************************************************************/
void destroyPomdpSolveParams( PomdpSolveParams param ) {
/* 
   Frees the memory for pointers in the params and the param structure
   itself.
*/
  if ( param->stat != NULL )
    destroySolutionStats( param->stat );

  if ( param->initial_policy != NULL )
    destroyAlphaList( param->initial_policy );

  free( param );

}  /* destroyPomdpSolveParams */
/**********************************************************************/
void parseMethodAliases( int argc, char **argv, int *mark_arg,
                         PomdpSolveParams param ) {

  /* zzz Add aliases here and see obsolete/old-parse-alias.c */

}  /* methodAliases */
/**********************************************************************/
void enforceSmallestPrecision( double *value, char *name ) {
/*
  Takes a value and makes sure it is not les than the smallest
  allowable precision the program uses.  It will give a message if
  it needs to be changed.  
*/
  char msg[MAX_MSG_LENGTH];

  if ( *value >= SMALLEST_PRECISION )
    return;
  
  *value = SMALLEST_PRECISION;
  
  sprintf( msg, 
           "The value for %s is below the smallest precision.\n\tSetting to %.3e." ,
           name, *value );

  Warning( msg )

}  /* enforceSmallestPrecision */
/**********************************************************************/
void parseCmdLineEpsilons( int argc, char **argv, int *mark_args,
                           PomdpSolveParams param) {

  /****************/
  /* Set precision of value function for each iteration when using the
     epsilon prune option. */
  getDoubleParam( argc, argv, 
                  CMD_ARG_PRUNE_EPSILON,
                  mark_args, 
                  &(param->prune_epsilon), 0.0, HUGE_VAL );

  enforceSmallestPrecision( &(param->prune_epsilon),
                            CMD_ARG_PRUNE_EPSILON );

  /* Set the main overall program precision. */
  getDoubleParam( argc, argv, 
                  CMD_ARG_EPSILON,
                  mark_args, 
                  &(param->epsilon), 0.0, HUGE_VAL );
    
  enforceSmallestPrecision( &(param->epsilon),
                            CMD_ARG_EPSILON );

  /* Set the precision to use in the LPs. */
  getDoubleParam( argc, argv, 
                  CMD_ARG_LP_EPSILON,
                  mark_args, 
                  &(param->lp_epsilon), 0.0, HUGE_VAL );
    
  enforceSmallestPrecision( &(param->lp_epsilon),
                            CMD_ARG_LP_EPSILON );

  /* Don't want to have the LPs be less precise than the rest of the
     program's operations. More precise is alright, since it can
     filter things out, but less precise makes little sense. */
  if ( param->lp_epsilon > param->epsilon ) {
    
    Warning( "LP epsilon must be no greater than general epsilon." );
    
    param->lp_epsilon = param->epsilon;

  } /* if lp_epsilon greater than general epsilon */

  /* Set the global LP precision
  LP_setPrecision( param->epsilon );

  /* This may need to change.  The relationship between these
     epsilons and the global epsilon might not be so
     straightforward. Also might need to add the
     enforceSmallestPrecision to be called on each of these if they
     are changed. */
  /* param->alpha_epsilon = param->epsilon; */
  param->vertex_epsilon = param->epsilon;
  param->double_equality_precision = param->epsilon;

  /* Because the impossible_obs_epsilon epsilon is related to
     probabilities in the model, an epsilon here near or greater than
     1 will mean that no observations will look possible.  For this
     reason, we have to take this to be some minimal value using
     epsilon nd default. */
  param->impossible_obs_epsilon = Min( param->epsilon, 
                                       DEFAULT_IMPOSSIBLE_OBS_EPSILON );
  
  /* Adjust the weak bound criteria, if none specified.  This
     requires us to look for the weak bound delta *after* the
     -epsilon argument so that we can override this if need
     be.  Added 8/17/95. */ 
  param->weak_bound_delta = param->epsilon;

  /****************/
  /* Set the weak bound delta (note it will not be used if we are not
     using the weak bound delta). */
  if ( getDoubleParam( argc, argv, 
                       CMD_ARG_DELTA,
                       mark_args, 
                       &(param->weak_bound_delta), 0.0, HUGE_VAL ))
    if ( param->stop_criteria == stop_weak )
      Warning
        ( "Weak bound delta on command line ignored.");
  
  enforceSmallestPrecision( &(param->weak_bound_delta),
                            CMD_ARG_DELTA );

} /* parseCmdLineEpsilons */
/**********************************************************************/
void parseCmdLineViVariations( int argc, char **argv, int *mark_args,
                               PomdpSolveParams param) {
/*
  There are some options for ways in which the value iteration can
  proceed.  These options have a bunch of tunable parameters,
  which are all parsed using this routine. 

  In particular one option is to try to maintain a fixed solution size
  for each epoch.  Another option is to adjust the precision of the
  algorithm as epochs evolve, only making things more precise when
  there is some form of convergence achieved.
*/

  /****************/
  /* The parameter defining what value iteration variation should be
     used. */
  getStringParamValidate( argc, argv, 
                          CMD_ARG_VI_VARIATION,
                          mark_args, 
                          (int *) &(param->vi_variation),
                          vi_variation_type_str, 
                          MAX_VI_VARIATION_TYPES );

  /****************/
  /* When using the adjustable epsilon value iteration variation, this
     sets the epsilon for the first epoch of value iteration. Also
     used for the fixed size solution variation. */
  getDoubleParam( argc, argv,  
                  CMD_ARG_STARTING_EPSILON,
                  mark_args, 
                  &(param->starting_epsilon), 0.0, 0.0 );

  /****************/
  /* When using the adjustable epsilon value iteration variation, this
     sets the epsilon for the final or most precise value. Also used
     for the fixed size variation of value iteration. */
  getDoubleParam( argc, argv, 
                  CMD_ARG_ENDING_EPSILON,
                  mark_args, 
                  &(param->ending_epsilon), 0.0, 0.0 );

  /****************/
  /* When using the adjustable epsilon value iteration variation, this
     sets the increment for the epsilon if it need to adjust. */
  getDoubleParam( argc, argv, 
                  CMD_ARG_EPSILON_ADJUST_FACTOR,
                  mark_args, 
                  &(param->epsilon_adjust_factor), 0.0, 0.0 );

  /****************/
  /* When using the fixed solution size value iteration variation,
     this sets the maximum number of vectors to allow as the solution
     on any epoch. Only parse this if that variation of VI was
     chosen. */
  if ( param->vi_variation == FixedSolnSizeVi ) {
    
    getIntParam( argc, argv, 
                 CMD_ARG_SOLN_SIZE,
                 mark_args, 
                 &(param->max_soln_size), 0, 0 );
    
    /* Make sure max solution size is not negative. */
    if ( param->max_soln_size < 1 ) {
      
      Warning( "max_soln is not positive.  Using normal VI." );
      param->max_soln_size = 0;
      param->vi_variation = NormalVi;
    } /* if max solution size makes no sense. */
     
  } /* if param->vi_variation == FixedSolnSizeVi */

  /****************/
  /* When using the adjustable epsilon value iteration variation, this
     helps define the criteria for when the epsilon should be
     adjusted.  It looks at the sizes of the resulting sets over the
     last few epochs and determines if they are all of approaximately
     the same size. This parameter defines how many epochs in the past
     it should look. */
  getIntParam( argc, argv, 
               CMD_ARG_HISTORY_WINDOW_LENGTH,
               mark_args, 
               &(param->epoch_history_window_length), 0, 0 );

  /****************/
  /* When using the adjustable epsilon value iteration variation, this
     helps define the criteria for when the epsilon should be
     adjusted. It looks at the sizes of the resulting sets over the
     last few epochs and determines if they are all of approximately
     the same size. This parameter defines what it means for the
     previous solution sizes to be approximately the same size.  They
     must all be within this many vectors of each other. */
  getIntParam( argc, argv,  
               CMD_ARG_HISTORY_WINDO_DELTA,
               mark_args, 
               &(param->epoch_history_window_delta), 0, 0 );

} /* parseCmdLineViVariations */
/**********************************************************************/
PomdpSolveParams parseCmdLineParams( int argc, 
                                     char **argv, 
                                     int *mark_args  ) {
/*
  Parses the command line for all the POMDP solution parameters.
*/
  PomdpSolveParams param;
  char str[MAX_CMD_ARG_LENGTH];
  int arg_num, i;

   /* We first create a parameter structure with the default values.
      Then  we'll go and read the command line args and over-ride
      whatever needs to be. */
   param = newPomdpSolveParams();

   /****************/
   /* Set brief reporting mode if specified.  Good for generating a 
      very succinct report when running a lot of experiments. */
   if ( getFlagParam( argc, argv, 
                      CMD_ARG_SUCCINCT,
                      mark_args ))
     param->succinct = TRUE;

   /****************/
   /* Whether or not to print out the statistic summary after the
      execution is completed. */
   if ( getFlagParam( argc, argv, 
                      CMD_ARG_STAT_SUMMARY,
                      mark_args ))
     param->stat_summary = TRUE;

   /****************/
   /* Set if we want to redirect everything to a file.  Note that we
      must do this early on and actually open the file here because we
      might shortly get error messages that will need to be printed
      out to the file. */
   if( getStringParam( argc, argv, 
                       CMD_ARG_REPORT_FILE,
                       mark_args, param->report_filename )) {
     
     if (( param->report_file 
          = fopen( param->report_filename , "w")) == NULL) {
       param->report_file = stdout;
       fprintf( stderr, 
                "** Error: Cannot write to output file %s.\n",
                param->report_filename );
       fprintf( stderr, 
                "\tUsing stdout instead.\n" );
     }  /* if can't open report file */

     /* If they desire to put all the output into a specific file,
        then we will also output all stderr messages here as well. */
     stderr = param->report_file;

   }  /* If report file specified */
   
   /* If no report file specified, so use the default 'stdout' */

   /****************/
   /* Set the name of the POMDP file to use. */
   if( !getStringParam( argc, argv, 
                        CMD_ARG_POMDP_FILE,
                        mark_args, param->param_filename ))

     /* If not found, indicate that there is no name, since we need to
        know this when we try to solve the problem.  We can't do
        anything without a POMDP file. */
     param->param_filename[0] = NULL_CHAR;

   /****************/
   /* See if we want to over-ride the discount factor.  It is
      usually specified in the file, but we want to allow it
      to be over-ridden on the command line. */
   getDoubleParam( argc, argv, 
                    CMD_ARG_DISCOUNT,
                    mark_args, 
                    &(param->override_discount), 0.0, 1.0 );

   /****************/
   /* Set which algorithm to use */
   arg_num = getStringParamValidate( argc, argv, 
                                     CMD_ARG_METHOD,
                                     mark_args, 
                                     (int *) &(param->method),
                                     method_str, MAX_NUM_METHODS );

   /* If there was a problem with getting the method, try some useful
      aliases which we allow. If no aliases are specified either, then
      this will return the default method. */
   if ( arg_num == 0 )
     parseMethodAliases( argc, argv, mark_args, param );

   /****************/
   /* Set which value iteration stopping criteria to use. */
   getStringParamValidate( argc, argv, 
                           CMD_ARG_STOP_CRITERIA,
                           mark_args, 
                           (int *) &(param->stop_criteria),
                           stop_criteria_str, 
                           MAX_NUM_STOP_CRITERIA );

   /****************/
   /* Set the horizon (if any) */
   getIntParam( argc, argv, 
                CMD_ARG_HORIZON,
                mark_args, 
                &(param->horizon), 0, 0 );

   /****************/
   /* See if we should set a timer. */
   getIntParam( argc, argv, 
                CMD_ARG_MAX_SECS,
                mark_args, 
                &(param->max_secs), 0, 0 );

   /****************/
   /* Set output file prefix to be used. */
   if( !getStringParam( argc, argv, 
                        CMD_ARG_OUTPUT,
                        mark_args, param->prefix_str ))
     /* If none specified, then use the default output file prefix */
     strcpy( param->prefix_str, DEFAULT_PREFIX );

   /****************/
   /* Set maximum amount of virtual memory it can use.  If the program
      attempts to use more than its limits, a SIGSEGV will be
      generated. */
   getIntParam( argc, argv, 
                CMD_ARG_MEMORY_LIMIT,
                mark_args, 
                &(param->memory_limit), 0, 0 );
   
   /****************/
   /* Set the name for the intitial policy file */
   if( !getStringParam( argc, argv, 
                        CMD_ARG_TERMINAL_VALUES,
                        mark_args, 
                        param->initial_policy_filename ))
     /* Else set the name to be an empty string.  This will indicate
        that we do not want to use an initial policy file. */
     param->initial_policy_filename[0] = NULL_CHAR;

   /****************/
   /* Save every set of vectors (after each epoch) */
   if ( getFlagParam( argc, argv, 
                      CMD_ARG_SAVE_ALL,
                      mark_args ))
     param->save_all = TRUE;

   /****************/
   /* Set how to purge projection sets. */
   getStringParamValidate( argc, argv, 
                           CMD_ARG_PROJ_PURGE,
                           mark_args, 
                           (int *) &(param->proj_purge),
                           purge_option_str, 
                           MAX_NUM_PURGE_OPTIONS );
   
   /****************/
   /* Set how to purge Q sets. */
   getStringParamValidate( argc, argv, 
                           CMD_ARG_Q_PURGE,
                           mark_args, 
                           (int *) &(param->q_purge_option),
                           purge_option_str, 
                           MAX_NUM_PURGE_OPTIONS );
   
   /****************/
   /* Set whether to use domination checks or not. */
   getStringParamValidate( argc, argv, 
                           CMD_ARG_DOM_CHECK,
                           mark_args, &(param->domination_check),
                           boolean_str, 
                           2 );
   
   /****************/
   /* Set whether to use witness poits or not. */
   getStringParamValidate( argc, argv, 
                           CMD_ARG_WITNESS_POINTS,
                           mark_args, &(param->use_witness_points),
                           boolean_str, 
                           2 );
   
   /****************/
   /* Set the number of random points to use in initialization of
      parsimonious sets for algorithms. */
   getIntParam( argc, argv, 
                CMD_ARG_ALG_INIT_RAND,
                mark_args, 
                &(param->alg_init_rand_points), 0, 0 );

   /****************/
   /* Set the number of random points to use in initialization of
      parsimonious sets for prune() routine. */
   getIntParam( argc, argv, 
                CMD_ARG_PRUNE_INIT_RAND,
                mark_args, 
                &(param->prune_init_rand_points), 0, 0 );

   /****************************************/
   /****  Algorithm specific section  ******/
   /****************************************/

   /* We must parse the '-method' argument first so we know which
      algorithm to send the parameters.  We could send them to all of
      them, but it seems wrong to let someone specify options which
      don't apply for an algorithm.  This way they will know if some
      options would be ignored. */

   /****************/
   /* See which version of incremental pruning to use */
   arg_num = getStringParamValidate( argc, argv, 
                                     CMD_ARG_INC_PRUNE_TYPE,
                                     mark_args, 
                                     (int *) &(param->ip_type),
                                     inc_prune_type_str, 
                                     MAX_INC_PRUNE_TYPES );

  /****************/
  /* Set how to purge the enumerated sets. */
  getStringParamValidate( argc, argv, 
                          CMD_ARG_ENUM_PURGE_OPTION,
                          mark_args, 
                          (int *) &(param->enum_purge_option),
                          purge_option_str, 
                          MAX_NUM_PURGE_OPTIONS );

   /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      The setting of the epsilons may need to change once I know 
      exactly how they all relate to each other.  For now, we set one
      global epsilon value and all the others are set to the same
      value. 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
   parseCmdLineEpsilons( argc, argv, mark_args, param );

   /* There are some options for ways in which the value iteration can
      proceed.  These options have a bunch of tunable parameters,
      which are all parsed using this routine. */
   parseCmdLineViVariations( argc, argv, mark_args, param );

  return ( param );

}  /* parseCmdLineParams */
/**********************************************************************/
void showPomdpSolveParams( PomdpSolveParams param ) {
/*
  Shows all the current settings and parameters for this run
  of the program.
*/
   int i;

   /* Don't show anything if we are in "succinct" mode. */
   if ( param->succinct == TRUE )
     return;

   fprintf( param->report_file, "       _________________\n");
   fprintf( param->report_file, "     <<<  POMDP-SOLVE  >>>\n");
   fprintf( param->report_file, "  ---------------------------\n");
   fprintf( param->report_file, "           PID=%d:\n", getPid() );

   /*******************************/
   /* General parameters          */
   /*******************************/
   fprintf( param->report_file, "General parameters:\n" );

   fprintf( param->report_file, "\tRandSeed = " );
   displayRandomSeed( param->report_file );
   fprintf( param->report_file, "\n" );

   /*******************************/
   /* Resource limits             */
   /*******************************/
   if (( param->max_secs > 0 ) || ( param->memory_limit > 0 ))
     fprintf( param->report_file, "Resource limit options:\n" );
   
   if ( param->max_secs > 0 )
     fprintf( param->report_file,"\tTime limit set to %d secs.\n", 
              param->max_secs );
   
   if ( param->memory_limit > 0 )
     fprintf( param->report_file,"\tMemory limit set to %d bytes.\n", 
              param->memory_limit );
   
   /*******************************/
   /* Value iteration parameters. */
   /*******************************/
   fprintf( param->report_file, "Value iteration parameters:\n" );

   fprintf( param->report_file, "\tPOMDP file = %s\n", 
            param->param_filename );

   if ( param->initial_policy_filename[0] != NULL_CHAR )
     fprintf( param->report_file, "\tInitial values = %s\n",
              param->initial_policy_filename );
   else
     fprintf( param->report_file, "\tInitial values = default\n" );

   if ( param->horizon < 1 )
     fprintf( param->report_file,"\tHorizon = Infinity\n");
   else
     fprintf( param->report_file,"\tHorizon = %d.\n", param->horizon );

   /* See if we over-rode the file's discount factor and if so,
      display what it is here. */
   if ( param->override_discount >= 0.0 ) 
     fprintf( param->report_file, "\tDiscount override = %e\n", 
              gDiscount );

   fprintf( param->report_file, "\tStopping criteria = %s", 
            stop_criteria_str[param->stop_criteria]);
   if ( param->stop_criteria == stop_weak )
     fprintf( param->report_file, " (delta = %e)\n", 
              param->weak_bound_delta );
   else
     fprintf( param->report_file, "\n" );

   /****************************************/
   /****  VI variation specific section  ***/
   /****************************************/

   fprintf( param->report_file, "\tVI Variation = %s\n", 
            vi_variation_type_str[param->vi_variation]);

   /* These are common to the adjustable epsilon and max solution size
      VI variations. */
   if (( param->vi_variation == AdjustableEpsilonVi ) 
       || ( param->vi_variation == FixedSolnSizeVi )) {
     
     fprintf( param->report_file, "\t\tStarting Epsilon = %e\n", 
              param->starting_epsilon );
     fprintf( param->report_file, "\t\tEnding Epsilon = %e\n", 
              param->ending_epsilon );
     fprintf( param->report_file, "\t\tEpsilon increment = %e\n", 
              param->epsilon_adjust_factor );

   } /* if AdjustableEpsilonVi or FixedSolnSizeVi */

  if ( param->vi_variation == AdjustableEpsilonVi ) {

     fprintf( param->report_file, "\t\tHistory window = %d\n", 
              param->epoch_history_window_length );
     fprintf( param->report_file, "\t\tHistory delta = %d\n", 
              param->epoch_history_window_delta );

   } /* if variation == AdjustableEpsilonVi */
   
  /*******************************/
  if ( param->vi_variation == FixedSolnSizeVi ) {
    
    fprintf( param->report_file,"\t\tMax Solution size: %d vectors.\n", 
             param->max_soln_size );
    
  } /* if variation == FixedSolnSizeVi */
   
  /*******************************/
  /* Optimization options.       */
  /*******************************/
  fprintf( param->report_file, "Optimization parameters:\n" );

   fprintf( param->report_file, "\tDomination check = %s\n", 
            boolean_str[param->domination_check]);

   fprintf( param->report_file, "\tGeneral Epsilon = %e\n", 
            param->epsilon );

   fprintf( param->report_file, "\tLP Epsilon = %e\n", 
            param->lp_epsilon );

   fprintf( param->report_file, "\tProjection purging = %s\n", 
            purge_option_str[param->proj_purge] );

   /* Only show the pruning epsilon if we are doing a non-normal
      prune of the Projection sets. */
   if ( param->proj_purge  == purge_epsilon_prune )
     fprintf( param->report_file, "\t\tPrune Epsilon = %e\n", 
              param->prune_epsilon );

   /* These options only apply to algorithms which build up their
      value functions one at a time, so don't show them for other
      algorithms. */
   switch ( param->method ) {
   case TwoPass:
   case Witness:
   case IncrementalPruning:
     fprintf( param->report_file, "\tQ purge = %s\n", 
              purge_option_str[param->q_purge_option]);

     /* Only show the pruning epsilon if we are doing a non-normal
        prune of the Q-sets. */
     if ( param->q_purge_option == purge_epsilon_prune )
       
       fprintf( param->report_file, "\t\tPrune Epsilon = %e\n", 
                param->prune_epsilon );
     break;
   case Enumeration:
   case LinearSupport:
   default:
     break;
   }  /* switch gMethod */
   

   fprintf( param->report_file, "\tUse witness points = %s\n", 
            boolean_str[param->use_witness_points]);

   if ( param->alg_init_rand_points > 0 )
     fprintf( param->report_file,
              "\tInit with %d random points.\n", 
              param->alg_init_rand_points );
   
   if ( param->prune_init_rand_points > 0 )
     fprintf( param->report_file,
              "\tInit pruning with %d random points.\n", 
              param->prune_init_rand_points );
   
  /*******************************/
  /* Algorithm options.          */
  /*******************************/
   fprintf( param->report_file, "Algorithm parameters:\n" );

   fprintf( param->report_file, "\tMethod = %s\n", 
            method_str[param->method]);

   /****************************************/
   /****  Algorithm specific section  ******/
   /****************************************/

   /* This defers to the actual individual modules themselves for
      algorithm-specific options. */
   switch ( param->method ) {
   case Enumeration:
     fprintf( param->report_file, "\tEnumeration method settings:\n" );
     fprintf( param->report_file, "\t\tEnumeration purging: %s\n", 
              purge_option_str[param->enum_purge_option] );
     break;
   case Witness:
     break;

   case TwoPass:
     break;

   case LinearSupport:
     break;

   case IncrementalPruning:
     fprintf( param->report_file, 
              "\tIncremental Pruning method settings:\n" );
     fprintf( param->report_file, "\t\tIncPrune type = %s\n", 
              inc_prune_type_str[param->ip_type]);
     break;

   default:
     break;
   }  /* switch gMethod */


   fprintf( param->report_file, "Solutions files:\n" );

   fprintf( param->report_file, "\tValues -> %s\n",
            param->alpha_filename );
   
   fprintf( param->report_file, "\tPolicy graph -> %s\n",
            param->pg_filename );
   
   if ( param->save_all == TRUE )
     fprintf( param->report_file, "\tSaving every epoch.\n" );
   
   fprintf( param->report_file, "  ---------------------------\n");
  
}  /* showPomdpSolveParams */
/**********************************************************************/
void showUsageParams( FILE *file ) {
  
  /*******************************/
  /* General parameters          */
  /*******************************/
  fprintf( file, "General options:\n" );

  fprintf( file, "\t%s <filename>\n", CMD_ARG_REPORT_FILE );

  fprintf( file, "\t%s <N1:N2:N3>\n", CMD_ARG_RAND_SEED );

  fprintf( file, "\t%s\n", CMD_ARG_SUCCINCT );
  
  fprintf( file, "\t%s\n", CMD_ARG_STAT_SUMMARY );
  
  /*******************************/
  /* Resource limits             */
  /*******************************/
  fprintf( file, "Resource limit options:\n" );

  fprintf( file, "\t%s <bytes>\n", CMD_ARG_MEMORY_LIMIT );

  fprintf( file, "\t%s <int>\n", CMD_ARG_MAX_SECS );

  /*******************************/
  /* Value iteration parameters. */
  /*******************************/
  fprintf( file, "Value iteration options:\n" );

  fprintf( file, "\t%s <pomdp-file>\n", CMD_ARG_POMDP_FILE );

  fprintf( file, "\t%s <policy-file>\n", CMD_ARG_TERMINAL_VALUES );

  fprintf( file, "\t%s <int>\n", CMD_ARG_HORIZON );

  fprintf( file, "\t%s (0-1)\n", CMD_ARG_DISCOUNT );

  showUsageEnumType( file,
                     CMD_ARG_STOP_CRITERIA,
                     MAX_NUM_STOP_CRITERIA,
                     stop_criteria_str );

  fprintf( file, "\t%s (0-infty)\n", CMD_ARG_DELTA );

  fprintf( file, "\t%s <file-prefix>\n", CMD_ARG_OUTPUT );

  fprintf( file, "\t%s\n", CMD_ARG_SAVE_ALL );

  /* Value iteration variation parameters. */
  showUsageEnumType( file,
                     CMD_ARG_VI_VARIATION,
                     MAX_VI_VARIATION_TYPES,
                     vi_variation_type_str );

  fprintf( file, "\t\t%s (0-infty)\n", CMD_ARG_STARTING_EPSILON );
  fprintf( file, "\t\t%s (0-infty)\n", CMD_ARG_ENDING_EPSILON );
  fprintf( file, "\t\t%s (0-infty)\n", CMD_ARG_EPSILON_ADJUST_FACTOR  );

  fprintf( file, "\t\t%s <int>\n", CMD_ARG_SOLN_SIZE );
  fprintf( file, "\t\t%s <int>\n", CMD_ARG_HISTORY_WINDOW_LENGTH );
  fprintf( file, "\t\t%s <int>\n", CMD_ARG_HISTORY_WINDO_DELTA );

  /*******************************/
  /* Optimization options.       */
  /*******************************/
  fprintf( file, "Optimization options:\n" );

  showUsageEnumType( file,
                     CMD_ARG_DOM_CHECK,
                     2,
                     boolean_str );

  fprintf( file, "\t%s (0-infty)\n", CMD_ARG_PRUNE_EPSILON );
  fprintf( file, "\t%s (0-infty)\n", CMD_ARG_EPSILON );
  fprintf( file, "\t%s (0-infty)\n", CMD_ARG_LP_EPSILON );

  showUsageEnumType( file,
                     CMD_ARG_PROJ_PURGE,
                     MAX_NUM_PURGE_OPTIONS,
                     purge_option_str );

  showUsageEnumType( file,
                     CMD_ARG_Q_PURGE,
                     MAX_NUM_PURGE_OPTIONS,
                     purge_option_str );

  showUsageEnumType( file,
                     CMD_ARG_WITNESS_POINTS,
                     2,
                     boolean_str );
  
  fprintf( file, "\t%s <int>\n", CMD_ARG_ALG_INIT_RAND );
  fprintf( file, "\t%s <int>\n", CMD_ARG_PRUNE_INIT_RAND );

  /*******************************/
  /* Algorithm options.          */
  /*******************************/
  fprintf( file, "Algorithm options:\n" );

  showUsageEnumType( file,
                     CMD_ARG_METHOD,
                     MAX_NUM_METHODS,
                     method_str );

  /****************************************/
  /****  Algorithm specific section  ******/
  /****************************************/

  fprintf( file, "  Enumeration method options:\n" );
  
  showUsageEnumType( file,
                     CMD_ARG_ENUM_PURGE_OPTION,
                     MAX_NUM_PURGE_OPTIONS,
                     purge_option_str );

  fprintf( file, "  Incremental pruning method options:\n" );

  showUsageEnumType( file,
                     CMD_ARG_INC_PRUNE_TYPE,
                     MAX_INC_PRUNE_TYPES,
                     inc_prune_type_str );

}  /* showUsageParams */
/**********************************************************************/
