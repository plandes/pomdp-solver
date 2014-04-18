/* 
  File: policy-graph.h
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
*/
#ifndef POLICY_GRAPH
#define POLICY_GRAPH

/**********************************************************************/
/********************       CONSTANTS       ***************************/
/**********************************************************************/

/**********************************************************************/
/********************   EXTERNAL VARIABLES   **************************/
/**********************************************************************/

/**********************************************************************/
/********************   EXTERNAL FUNCTIONS    *************************/
/**********************************************************************/

/* Displays the policy graph to the file handle specified.

  The policy graph will be output with the format of one line per node
  in the policy graph:

  ID  ACTION    OBS1  OBS2 OBS3 ... OBSN

  where ID is the id of the alpha vector in the current set, ACTION is
  the action for this vector, and OBS1 through OBSN are the id's of
  the vectors in the previous epoch's alpha vector set (one for each
  observation).  */
extern void displayPolicyGraph( FILE *file, AlphaList list );

/* Displays the policy graph of a set of vectors to the filename
  specified.  */
extern void writePolicyGraph( AlphaList list, char *filename );

/* Displays the policy graph for an alpha vector set to stdout.  */
extern void showPolicyGraph( AlphaList list );
   
#endif
