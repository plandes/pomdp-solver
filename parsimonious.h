/* 
  File: parsimonious.h
  Author: A. R. Cassandra
  July, 1998

  *****
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
*/

#ifndef PARSIMONIOUS_H
#define PARSIMONIOUS_H

#include "params.h"

/**********************************************************************/
/********************       CONSTANTS       ***************************/
/**********************************************************************/

/**********************************************************************/
/********************   DEFAULT VALUES       **************************/
/**********************************************************************/

/**********************************************************************/
/********************   EXTERNAL VARIABLES   **************************/
/**********************************************************************/

/**********************************************************************/
/********************   EXTERNAL FUNCTIONS    *************************/
/**********************************************************************/

/* Removes all vectors from the list that can be determined to have a
  non-empty region using only the simple component-wise domination
  check with some other vector in the list.  */
extern int dominationCheck( AlphaList orig_list );

/* Will generate 'num_points' random belief points and mark the
   vectors at these points to the list.  */
extern void markBestAtRandomPoints( AlphaList list, 
                                    int num_points, 
                                    int save_witness_points,
                                    double epsilon );

/* Sets the 'mark' field of each vector that dominates at some belief
   simplex vertex.  A difference between this routine and
   initWithSimplexCornersQ() is that this vector does not have to
   construct the vectors for the points, it simply picks them out of
   the list provided.  This is used in the prune algorithm to
   initialize the list.  Will loop through all the belief simplex
   vertices and find the best vectors at these points from the list.  */
extern void markBestAtSimplexVertices( AlphaList list, 
                                int save_witness_points, 
                                double epsilon );

/* There are a number of simple checks that can be made to determine
  if the region will be non-empty.  First, if the 'alpha' vector is
  already in the list.  Second, if something in the list
  component-wise dominates the vector.  This routine does these checks
  and returns TRUE if region can easily be determined to be empty and
  FALSE if the simple check reveals nothing.  */
extern int isEmptyRegionSimpleCheck( AlphaList list, 
                                     double *alpha,
                                     double epsilon,
                                     int domination_check );
 
/* Will use linear programming to prune the list to a unique
   parsimonious represenation.  If the save_points flag is set, then
   each vector in the resulting set (with a non-empty) region will
   have the witness point used to verify its non-empty region saved in
   the node containing the vector. Returns the number of nodes pruned.

   This uses the scheme propose by Lark and White, which is mentioned
   in White's 1991 Operations Research POMDP survey article.  */
extern int prune( AlphaList orig_list, 
                  PurgeOption purge_option,
                  PomdpSolveParams param );
     
/* Removes vectors from the list according to the purging option sent
  it.  It can do anything from nothing, to simple domination checks,
  to full blown 'pruning'.  */
extern void purgeAlphaList( AlphaList list, 
                            PurgeOption purge_option,
                            PomdpSolveParams param );
   
/* Runs the purgeAlphaList() routine on all the projection sets using
   the purging option for projections as set on the command line (or
   with default.)  */
extern void purgeProjections( AlphaList **projection, 
                              PomdpSolveParams param );

#endif
