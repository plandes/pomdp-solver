/* 
  File: witness.h
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
#ifndef WITNESS_H
#define WITNESS_H

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

extern void initWitness( );

extern void cleanUpWitness( );

/* The main witness algorithm routine for finding the Q-function
  represention for value iteration with POMDPs.  */
extern AlphaList improveWitness( AlphaList *projection,
                                 PomdpSolveParams param );
      
#endif
