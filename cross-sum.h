/*
  File: cross-sum.h
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
  ***
*/
#ifndef CROSS_SUM_H
#define CROSS_SUM_H

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

/* Takes the cross sum of two sets of vectors and returns the
  resulting set. If either A or B is null, then NULL is returned.  If
  either list is empty, then an empty list is returned. The
  save_obs_sources argument deterines whther we do the bookkeeping
  required to develop a policy graph or not.  */
extern AlphaList crossSum( AlphaList A, 
                           AlphaList B, 
                           int save_obs_sources );
 
#endif
