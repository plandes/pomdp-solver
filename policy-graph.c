/*
  File: policy-graph.c
  Author: Anthony R. Cassandra
  July, 1998

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

  Routines that deal with policy graphs.
*/

#define POLICY_GRAPH_C

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <mdp.h>

#include "global.h"
#include "pomdp.h"
#include "alpha.h"
#include "policy-graph.h"

/**********************************************************************/
void displayPolicyGraph( FILE *file, AlphaList list ) {
/*
  Displays the policy graph to the file handle specified.

  The policy graph will be output with the format of one line per node
  in the policy graph:

  ID  ACTION    OBS1  OBS2 OBS3 ... OBSN

  where ID is the id of the alpha vector in the current set, ACTION is
  the action for this vector, and OBS1 through OBSN are the id's of
  the vectors in the previous epoch's alpha vector set (one for each
  observation). 
*/
  int z;

  Assert( file != NULL && list != NULL,
          "Bad (NULL) parameter(s)." );

  list = list->head;

  while( list != NULL ) {

    fprintf( file, "%d %d    ", list->id, list->action );

    if ( list->obs_source != NULL )
      for ( z = 0; z < gNumObservations; z++ ) {

        if ( list->obs_source[z] != NULL )
          fprintf( file, " %4d", list->obs_source[z]->id );
        else
          /* We put an 'X' when that observation is impossible. */
          fprintf( file, "    X" );
      }  /* for z */          

    else
      fprintf( file, "[No information available]" );

    fprintf( file, "\n" );
    
    list = list->next;
  }  /* while */

}  /* displayPolicyGraph */
/**********************************************************************/
void writePolicyGraph( AlphaList list, char *filename ) {
/*
  Displays the policy graph of a set of vectors to the filename
  specified. 
*/
   FILE *file;

   if ((file = fopen(filename , "w")) == NULL) {
     fprintf( stderr, 
             "** Error: The policy graph file: %s cannot be opened.\n",
             filename);
      return;
   }

   displayPolicyGraph( file, list );
   
   fclose( file );

}  /* writePolicyGraph */
/**********************************************************************/
void showPolicyGraph( AlphaList list ) {
/*
  Displays the policy graph for an alpha vector set to stdout.
*/
  displayPolicyGraph( stdout, list );
}  /* showPolicyGraph */
/**********************************************************************/
