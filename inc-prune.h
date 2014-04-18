/* 
  File: inc-prune.h
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
#ifndef INC_PRUNE_H
#define INC_PRUNE_H

#include "params.h"

/**********************************************************************/
/********************       CONSTANTS       ***************************/
/**********************************************************************/

/* We will use the 'length' field of an AlphaList as a counter.  To
   make this explicit, this macros is used. */
#define COUNT(X)    ((X)->length)

/* Use the 'mark' field of the header node of a list to indicate
   whether or not it will need to be destroyed later. The use of this
   is made explicit with this macro. */
#define SHOULD_DESTROY(X) ((X)->mark)

/**********************************************************************/
/********************   DEFAULT VALUES       **************************/
/**********************************************************************/

/**********************************************************************/
/********************   EXTERNAL VARIABLES   **************************/
/**********************************************************************/

/**********************************************************************/
/********************   EXTERNAL FUNCTIONS    *************************/
/**********************************************************************/

extern void initIncPrune( );

extern void cleanUpIncPrune(  );
  
/* The main incremental pruning algorithm routine for finding the
  Q-function represention for value iteration with POMDPs.  */
extern AlphaList improveIncPrune( AlphaList *projection,
                                  PomdpSolveParams param );
  
#endif
