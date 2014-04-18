/* 
  File: global.c 
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <mdp.h>

#include "global.h"

/* Nice to have pretty names for the boolean values 0 and 1. */
char *boolean_str[] = BOOLEAN_STRINGS;

/* Strings for the various stopping criteria */
char *verbose_mode_str[] = VERBOSE_MODE_STRINGS;

static FILE *gStdErrFile = stderr;

/* The name of the executable of this program. */
char gExecutableName[80];

/* There are various ways to turn verboseness on and off.  Each elemnt
   of the array defines whether one of these is on or off.  The
   mnemonics in the header file show which ones are which. */
int gVerbose[NUM_VERBOSE_MODES];

/**********************************************************************/
/* Temporary variables usefule for scratch work. */
/**********************************************************************/

/* There are times when we need an array of doubles for temporary
   workspace.  These vectors will be gNumStates in length. */
double *gTempValue;
double *gTempBelief;
double *gTempAlpha;

/**********************************************************************/
void initGlobal(  ) {
/*
  Sets up and allocates variables that are used globally across
  modules in the program. Currently just allocates a bunch of scratch
  memory areas.
*/
  gTempBelief = (double *) malloc( gNumStates * sizeof( double ));
  gTempAlpha = (double *) malloc( gNumStates * sizeof( double ));
  gTempValue = (double *) malloc( gNumStates * sizeof( double ));
  
}  /* initGlobal */
/**********************************************************************/
void cleanUpGlobal(  ) {
/*
  Cleans up after problem is solved to free any resources and reset
  anything that the initGlobal() routine did.
*/

  free( gTempBelief );
  free( gTempAlpha );
  free( gTempValue );

}  /* cleanUpGlobal */
/**********************************************************************/
int getPid(  ) {
/* 
   Just a wrapper to the UN*X getpid() function to isolate it in case
   this gets ported to another platform.  Note that for POSIX, the
   'pid_t' type returned by getpid() is an 'int'.
*/
  return( (int) getpid() );
}  /* getPid */
/**********************************************************************/
void removeFile( char *filename ) {
/* 
   Just a wrapper to the UN*X unlink() function to isolate it in case
   this gets ported to another platform.  
*/

  unlink( filename );

}  /* removeFile */
/**********************************************************************/
