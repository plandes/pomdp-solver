/* pomdp.h

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
*/
#ifndef POMDP_H
#define POMDP_H

#include "global.h"
#include <mdp.h>

/**********************************************************************/
/********************       CONSTANTS       ***************************/
/**********************************************************************/

/**********************************************************************/
/********************   EXTERNAL VARIABLES   **************************/
/**********************************************************************/

extern int **gObservationPossible;
extern int *gNumPossibleObservations;

/**********************************************************************/
/********************   EXTERNAL FUNCTIONS    *************************/
/**********************************************************************/

/* Often we would like to do some max or min procedure and require
   initialization to the most extreme value.  Since the extreme value
   depends on whether or not we are using rewards or costs, we have
   encapsulated this in this routine.  */
extern double worstPossibleValue();

/* Often we would like to do some max or min procedure and require
   initialization to the most extreme value.  Since the extreme value
   depends on whether or not we are using rewards or costs, we have
   encapsulated this in this routine.  */
extern double bestPossibleValue();

/* Often we would like to do some max or min procedure and require
   coparing a new value to a current value.  Since the test for which
   is better depends on whether rewards or costs are being used, we
   have encapsulated this in this routine.  We also want to account
   for the precision of the current run (i.e.,
   gDoubleEqualityPrecision.)  */
extern int isBetterValue( double new_value, 
                          double current, 
                          double epsilon );

/* Does the necessary things to read in and set-up a POMDP file.  Also
 precomputes which observations are possible and which are not.  */
extern void initializePomdp( char *filename, 
                             double obs_possible_epsilon );
    
/* Deallocates the POMDP read in by initializePomdp().  */
extern void cleanUpPomdp(  );
 
#endif
