/*  
  File: pomdp.c
  Author: Anthony R. Cassandra
  July, 1998

  *****
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

  This file contains code for reading in a pomdp file and setting
  the global variables for the problem for use in all the other files.
  It also has routines for operations on belief states.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mdp.h>

#include "global.h"
#include "pomdp.h"

int **gObservationPossible;
int *gNumPossibleObservations;

/**********************************************************************/
double worstPossibleValue() {
/*
   Often we would like to do some max or min procedure and require
   initialization to the most extreme value.  Since the extreme value
   depends on whether or not we are using rewards or costs, we have
   encapsulated this in this routine.  */

   if( gValueType == REWARD_value_type )
     return( -1.0 * HUGE_VAL );
   else
     return( HUGE_VAL );

}  /* worstPossibleValue */
/**********************************************************************/
double bestPossibleValue() {
/*
   Often we would like to do some max or min procedure and require
   initialization to the most extreme value.  Since the extreme value
   depends on whether or not we are using rewards or costs, we have
   encapsulated this in this routine.  */

   if( gValueType == REWARD_value_type )
     return( HUGE_VAL );
   else
     return( -1.0 * HUGE_VAL );

}  /* bestPossibleValue */
/**********************************************************************/
int isBetterValue( double new_value, double current, double epsilon ) {
/* 
   Often we would like to do some max or min procedure and require
   coparing a new value to a current value.  Since the test for which
   is better depends on whether rewards or costs are being used, we
   have encapsulated this in this routine.  We also want to account
   for the precision of the current run (i.e.,
   gDoubleEqualityPrecision.) 
*/

  if( gValueType == REWARD_value_type )
    return( LessThan( current, new_value, epsilon ));
  else
    return( LessThan( new_value, current, epsilon ));

}  /* isBetterValue */
/**********************************************************************/
void setPossibleObservations( double epsilon ) {
/*
  Sets the global arrays to precomputed values to determine whether or
  not each observation is possible for a given action.  Also stores
  how many observations are possible for each action.
*/
  int a, z, j, cur_state;
  int all_zero_prob_obs;

  for ( a = 0; a < gNumActions; a++ ) {

    for ( z = 0; z < gNumObservations; z++ ) {
      
      /* We want to check for the case where an observation is
         impossible.  */

      all_zero_prob_obs = TRUE;
      for ( cur_state = 0; cur_state < gNumStates; cur_state++)
        for ( j = P[a]->row_start[cur_state]; 
              j < P[a]->row_start[cur_state] 
                + P[a]->row_length[cur_state];
              j++ ) 
          if ( ! Equal( getEntryMatrix( R[a], P[a]->col[j], z ),
                        0.0, epsilon )) {
            all_zero_prob_obs = FALSE;
      
            /* Yeah, it's a 'goto'; just so I can say I used one. */
            goto END_LOOP;
          }
      
    END_LOOP:
      
      if ( all_zero_prob_obs )
        gObservationPossible[a][z] = FALSE;
      
      else  {
        gObservationPossible[a][z] = TRUE;
        gNumPossibleObservations[a]++;
      }  /* if observation is possible */
      
    } /* for z */
  
  } /* for a */

  /* A little sanity check. */
  for ( a = 0; a < gNumActions; a++ )
    Assert( gNumPossibleObservations[a] > 0,
            "Bad POMDP. No observations possible for some action." );

}  /* setPossibleObservations */
/**********************************************************************/
void initializePomdp( char *filename, 
                      double obs_possible_epsilon ) {
/*
  Does the necessary things to read in and set-up a POMDP file.
  Also precomputes which observations are possible and which are not.
*/
  int a;
  char msg[MAX_MSG_LENGTH];
  
  if (( filename == NULL ) || ( filename[0] == NULL_CHAR )) {
    sprintf( msg,
             "No parameter file specified (Use '%s' for options.)",
             CMD_ARG_HELP_SHORT );
    Abort( msg );
  }
  
  if ( ! readMDP( filename )) {
    sprintf( msg, "Could not successfully parse file: %s.\n",
             filename );
    Abort( msg );
  } /* if problem parsing POMDP file. */
  
  /* The code no longer supports anything but rewards and I should
     make this explicit.  */
  if ( gValueType == COST_value_type ) {
    sprintf( msg, "%s\n\t%s",
             "This program no longer supports problems with costs.",
             "Multiply costs by '-1' and use rewards." );
    Abort( msg );

  } /* if rewards are specified. */

  
  if ( gProblemType != POMDP_problem_type ){
    sprintf( msg,
             "Parameter file is not a POMDP specification." );
    Abort( msg );
  }

  /* We'll use this stuff if the setPossibleObservations() routine is
     called. */ 
  gObservationPossible 
    = (int **) malloc( gNumActions 
                       * sizeof( *gObservationPossible ));
  for ( a = 0; a < gNumActions; a++ )
    gObservationPossible[a]
      = (int *) calloc( gNumObservations, 
                        sizeof( **gObservationPossible ));
  
  gNumPossibleObservations
    = (int *) calloc( gNumActions,
                      sizeof( *gNumPossibleObservations ));

  setPossibleObservations( obs_possible_epsilon );

}  /* initializePomdp */
/**********************************************************************/
void cleanUpPomdp(  ) {
/*
  Deallocates the POMDP read in by initializePomdp().
*/
  int a;

  for ( a = 0; a < gNumActions; a++ )
    free( gObservationPossible[a] );
  free( gObservationPossible );
  
  free( gNumPossibleObservations );
  
  deallocateMDP();

}  /* cleanUpPomdp */
/**********************************************************************/
