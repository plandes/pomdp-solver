/* 
  File: timing.h
  Author: Anthony R. Cassandra
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
#ifndef TIMING_H
#define TIMING_H

/*******************************************************************/
/**************       EXTERNAL VARIABLES            ****************/
/*******************************************************************/

/*******************************************************************/
/**************       EXTERNAL FUNCTIONS            ****************/
/*******************************************************************/

/* Get total CPU time in seconds breaking it down by user and system
   time. */
extern void getSecsDetail( double *user_time, 
                           double *system_time );

/* Get total CPU time in seconds including user and system time. */
extern double getSecs( );

/* Report the total secons time in a nicer hr, min sec format with a
   string to label what the time is for.  */
extern void reportTimes( FILE *file,
                         double tot_secs, 
                         char *str );
   
#endif
