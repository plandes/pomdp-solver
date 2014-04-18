/*
  File: auto-test.c
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

  Runs a series of tests and prints out places where it fails.  Note
  that an error in one routine might result in many routines seeming
  to have errors.

  It is not a completely exhaustive test, but exercises the basic
  functionality.  As bugs are found, test cases should be added to
  detect that they have been fixed.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Note that the ordering of the inclusion of these header files
   roughly reflects the dependencies between them, */

#include <mdp.h>

#include "global.h"
#include "timing.h"
#include "random.h"
#include "pomdp.h"
#include "alpha.h"
#include "params.h"
#include "stats.h"
#include "projection.h"
#include "lp-interface.h"
#include "region.h"
#include "parsimonious.h"
#include "neighbor.h"
#include "common.h"
#include "cross-sum.h"

#include "enumeration.h"
#include "linear-support.h"
#include "two-pass.h"
#include "witness.h"
#include "inc-prune.h"

#include "signal-handler.h"
#include "policy-graph.h"
#include "pomdp-solve.h"

#define error(MSG) { \
  printf( "** Error ** " ); \
  printf( MSG ); \
  printf( "\n" ); \
  gTotErrors++; \
}

/* Used if you want a more verbose output during testing. */
int gAutoTestVerbose = FALSE;

int gTotErrors = 0;

/********************* random.c test  *****************************/
/* How many slots to use in random number test. */
#define NUM_RANDOM_SLOTS               100

/* How many random numbers to generate in randmo number test. */
#define NUM_RANDOM_NUMBERS            100000

#define RAND_DOUBLE_MIN                -10.8
#define RAND_DOUBLE_MAX                +25.5
#define RAND_INT_MIN                   -47

#define LP_TEMP_FILENAME               "zzz.lp.zzz"
#define SOLN_TEMP_FILENAME             "zzz.out.zzz"

#define LP_EPSILON                     1e-9

/* When lp_solve writes the solution out to a file, it uses 5 decimal
   places.  We aren't too concerned with precision here, so we do not
   force lp_solve to print more precision, we just use a larger
   epsilon when comparing the internal and external solutions. */
#define LP_SOLN_EPSILON                1e-3

/********************* lp-interface.c test  *********************/

#define NUM_RANDOM_LPS                 25
#define RAND_LP_MIN_ROWS                5
#define RAND_LP_MAX_ROWS               30
#define RAND_LP_MIN_COLS                5
#define RAND_LP_MAX_COLS               30
#define MIN_LP_COEF                     0.0
#define MAX_LP_COEF                    10.0

int lp_rows, lp_cols, lp_obj_sense;

double *lp_obj_coef;
double **lp_coef_matrix;
double *lp_rhs;
char *lp_sense, *lp_sense_out_char;
int lp_feasible;
int *lp_row_has_coef;

double lp_obj_solution;
double *lp_solution;
int lp_has_equality;

/*********************    region.c test      *********************/

#define REGION_TEST_STATES             15

/********************* parsimonious.c test  *********************/

#define PARSIMONIOUS_TEST_STATES       10
#define PARSIMONIOUS_RAND_ALPHAS       20

/* To encode values so checking can be done. Two must be consistent. */
#define NEIGHBOR_NUM_ROWS                   8
#define NEIGHBOR_ROW_BITS                   3

/********************* cross-sum.c test  *********************/

#define NUM_ENCODED_ALPHA_VECTORS          16
#define NUM_ENCODED_SHIFT_BITS              4

/********************* Algorithms test      *********************/

#define NUM_SOLVING_TESTS                   10

/* This defines which algorithms should be checked. */
static int AlgorithmTesting[MAX_NUM_METHODS] = {
  TRUE,    /* Enumeration */
  TRUE,    /* Twopass */
  FALSE,    /* LinearSupport */
  TRUE,     /* Witness */
  TRUE     /* IncrementalPruning */
  
};

static int NumAlgorithmVariations[MAX_NUM_METHODS] = {
  1,    /* Enumeration */
  1,    /* Twopass */
  1,    /* LinearSupport */
  1,    /* Witness */
  3     /* IncrementalPruning */
  
};

static char *SolveDescriptStr[NUM_SOLVING_TESTS] = {
  "machine maintenance (initial)",
  "machine maintenance",
  "aaai tiger",
  "cheese maze",
  "4x4",
  "part painting",
  "shuttle",
  "4x3",
  "network",
  "aircraft ident"
};

static char *SolveParamStr[NUM_SOLVING_TESTS] = {
  "solutions/ejs1.POMDP",
  "solutions/ejs1.POMDP",
  "solutions/tiger.95.POMDP",
  "solutions/cheese.95.POMDP",
  "solutions/4x4.95.POMDP",
  "solutions/hanks.95.POMDP",
  "solutions/shuttle.95.POMDP",
  "solutions/4x3.95.POMDP",
  "solutions/network.POMDP",
  "solutions/saci-s12-a6-z5.POMDP"
};

static char *SolveInitValStr[NUM_SOLVING_TESTS] = {
  "solutions/ejs1.initial.alpha",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  ""
};

static char *SolveSolnStr[NUM_SOLVING_TESTS] = {
  "solutions/ejs1.init.alpha",
  "solutions/ejs1.alpha",
  "solutions/tiger.alpha",
  "solutions/cheese.95.alpha",
  "solutions/4x4.95.alpha",
  "solutions/hanks.95.alpha",
  "solutions/shuttle.95.alpha",
  "solutions/4x3.95.alpha",
  "solutions/network.95.alpha",
  "solutions/saci-s12-a6-z5.95.alpha"
};

/* Because different algorithms take varying amounts of time on
   different problems, we individually say how long to run each for so
   that the test run does not take too long. */ 
static int SolveHorizon[NUM_SOLVING_TESTS][MAX_NUM_METHODS] = {

/*-----------------------------------------------------------
 enum       two-pass    lin-sup     witness     inc-prune 
-------------------------------------------------------------*/
  20,         20,         20,         20,         20,   /* ejs init */
  20,         20,         20,         20,         20,   /* ejs */
  15,         10,         15,         15,         15,   /* tiger */
   5,          5,          5,          5,          5,   /* cheese */
   5,          5,          5,          5,          5,   /* 4x4 */
   5,          5,          5,          5,          5,   /* hanks */
   5,          5,          5,          5,          5,   /* shuttle */
   5,          5,          5,          5,          5,   /* 4x3 */
   3,          3,          3,          3,          3,   /* network */
   3,          3,          3,          3,          2    /* saci */
};

/********************* Policy Graph test      *********************/

#define NUM_PG_TESTS                   4

static char *PgDescriptStr[NUM_PG_TESTS] = {
  "aaai tiger",
  "cheese maze",
  "4x4",
  "part painting"
};

static char *PgParamStr[NUM_PG_TESTS] = {
  "solutions/tiger.95.POMDP",
  "solutions/cheese.95.POMDP",
  "solutions/4x4.95.POMDP",
  "solutions/hanks.95.POMDP"
};

static char *PgInitValStr[NUM_PG_TESTS] = {
  "solutions/tiger.95.optimal.alpha",
  "solutions/cheese.95.optimal.alpha",
  "solutions/4x4.95.optimal.alpha",
  "solutions/hanks.95.optimal.alpha"
};

static char *PgSolnStr[NUM_PG_TESTS] = {
  "solutions/tiger.95.optimal.pg",
  "solutions/cheese.95.optimal.pg",
  "solutions/4x4.95.optimal.pg",
  "solutions/hanks.95.optimal.pg"
};

/**********************************************************************/
/**********    UTILITY ROUTINES         *******************************/
/**********************************************************************/

/**********************************************************************/
AlphaList newIdentityAlphaList(  ) {
/*
  This sets up a list that would look like an identity matrix when
  viewing alpha vector components as the cols and the vectors as
  the rows. 
*/ 
  int i, j;
  AlphaList list;
  double *alpha;
  
  list = newAlphaList();

  for ( i = 0; i < gNumStates; i++ ) {
    alpha = newAlpha();
    for ( j = 0; j < gNumStates; j++ )
      if ( i == j )
        alpha[j] = 1.0;
      else
        alpha[j] = 0.0;
    appendAlphaList( list, alpha, 0 );
  } /* for i */

  return ( list );

}  /* newIdentityAlphaList */
/**********************************************************************/
AlphaList createAlphaList( int n ) {
/*
  Creates an AlphaList with 'n' nodes.  The first node is all 1's,
  the second is all 2's and the third is all 3's, etc.
*/
  AlphaList list;
  double *alpha;
  int i, j;

  list = newAlphaList();
  
  for ( i = 0; i < n; i++ ) {
    alpha = newAlpha();
    for ( j = 0; j < gNumStates; j++ ) 
      alpha[j] = i;

    appendAlphaList( list, alpha, 0 );
    
  } /* for i */

  return ( list );

}  /* createThreeNodeList */
/**********************************************************************/
void createRandomLP( int sparse ) {
  int r, c;

  lp_has_equality = FALSE;

  /* Get a random size for the LP */
  lp_rows = getRandomInt( RAND_LP_MIN_ROWS, RAND_LP_MAX_ROWS );
  lp_cols = getRandomInt( RAND_LP_MIN_COLS, RAND_LP_MAX_COLS );

  /* Allocate space for the LP (dense). */
  lp_obj_coef = (double *) malloc( lp_cols * sizeof( double ));

  lp_coef_matrix = (double **) malloc( lp_rows * sizeof( *lp_coef_matrix ));
  for ( r = 0; r < lp_rows; r++ )
    lp_coef_matrix[r] = (double *) malloc( lp_cols * sizeof( double ));

  lp_rhs = (double *) malloc( lp_rows * sizeof( double ));
  lp_sense = (char *) malloc( lp_rows * sizeof( char ));
  lp_sense_out_char = (char *) malloc( lp_rows * sizeof( char ));
  lp_row_has_coef = (int *) calloc( lp_rows, sizeof( int ));

  /* Set up the LP */
  if ( fran() < 0.5 )
    lp_obj_sense = MAXIMIZE;
  else
    lp_obj_sense = MINIMIZE;

  for ( c = 0; c < lp_cols; c++ ) 
    lp_obj_coef[c] = getRandomDouble( MIN_LP_COEF, MAX_LP_COEF );

  /* Set coef matrix (make about half the entries zero if sparse is
     specified.) A tricky part here is that we must ensure that each
     row has at least one non-zero coefficient. We do this after this
     doubly-nested loop. */
  for ( r = 0; r < lp_rows; r++ )
    for ( c = 0; c < lp_cols; c++ ) 
      if ( (! sparse) || (fran() < 0.5) ) {
        lp_coef_matrix[r][c] 
          = getRandomDouble( MIN_LP_COEF, MAX_LP_COEF );
        lp_row_has_coef[r] = TRUE;
      }
      else
        lp_coef_matrix[r][c] = 0.0;

  /* For any row that didn't have at least one non-zero row set, pick
     a random column and set it to a random value. */
  for ( r = 0; r < lp_rows; r++ )
    if ( ! lp_row_has_coef[r] )
      lp_coef_matrix[r][getRandomInt( 0, lp_cols-1 )]
        = getRandomDouble( MIN_LP_COEF, MAX_LP_COEF );

  for ( r = 0; r < lp_rows; r++ )
    lp_rhs[r] = getRandomDouble( MIN_LP_COEF, MAX_LP_COEF );

  for ( r = 0; r < lp_rows; r++ )
    if ( ! lp_has_equality
         && ( fran() < 0.05 )) {
      lp_sense[r] = 'E';
      lp_sense_out_char[r] = '=';
      lp_has_equality = 1;
    }
  
  /* We'll deterministically use <= for MAX and >= for MIN problems. */ 
    else if ( lp_obj_sense == MINIMIZE ) {
      lp_sense[r] = 'G';
      lp_sense_out_char[r] = '>';
    }
    else {
      lp_sense[r] = 'L';
      lp_sense_out_char[r] = '<';
    }

  lp_solution = (double *) malloc( lp_cols * sizeof( double ));

}  /* createRandomLP */
/**********************************************************************/
void freeRandomLP(  ) {
  int r;

  lp_rows = lp_cols = 0;

  free( lp_obj_coef );

  for ( r = 0; r < lp_rows; r++ )
    free( lp_coef_matrix[r] );
  free( lp_coef_matrix );

  free ( lp_rhs );
  free( lp_sense );
  free( lp_sense_out_char );
  free ( lp_row_has_coef );

  free( lp_solution );
}  /* freeRandomLP */
/**********************************************************************/
LP setFromRandomLP( int sparse ) {
  LP lp;
  int r, c;
  int non_zeroes = 0;
  int index = 0;

  if ( ! sparse )
    non_zeroes = lp_rows * lp_cols;
  else
    for ( r = 0; r < lp_rows; r++ )
      for ( c = 0; c < lp_cols; c++ ) 
        if ( ! Equal( lp_coef_matrix[r][c], 0.0, LP_EPSILON ))
          non_zeroes++;
  
  lp = LP_newLP( lp_rows, lp_cols, non_zeroes );

  lp->objsen = lp_obj_sense;
  for ( c = 0; c < lp_cols; c++ ) {
    lp->obj[c] = lp_obj_coef[c];
    lp->lowbnd[c] = 0.0;
    lp->upbnd[c] = INFBOUND;
  } /* for i */

  for ( c = 0; c < lp_cols; c++ ) {
    lp->matbeg[c] = index;
    for ( r = 0; r < lp_rows; r++ ) {
      if ( sparse 
           && Equal( lp_coef_matrix[r][c], 0.0, LP_EPSILON )) 
        continue;
      
      lp->matval[index] = lp_coef_matrix[r][c];
      lp->matind[index++] = r;
    } /* for r */
    lp->matcnt[c] = index - lp->matbeg[c];
  } /* for c */

  if ( index != lp->matspace )
    error( "** Error ** LP allocated non-zero diff from actual." );

  for ( r = 0; r < lp_rows; r++ ) {
    lp->rhs[r] = lp_rhs[r];
    lp->sense[r] = lp_sense[r];
  } /* for r */
      
  return ( lp );

}  /* setFromRandomLP */
/**********************************************************************/
void writeExternalCplex(  ) {
  Abort( "writeExternalCplex() not yet implemented." );
} /* writeExternalCplex */
/**********************************************************************/
void writeExternalLpSolve(  ) {
  FILE *file;
  int r,c;

  if ((file = fopen( LP_TEMP_FILENAME, "w")) == NULL) 
    Abort( "Cannot open LP temp file for writing." );

  if ( lp_obj_sense == MAXIMIZE )
    fprintf( file, "max: " );
  else
    fprintf( file, "min: " );

  for ( c = 0; c < lp_cols; c++ ) 
    fprintf( file, "%+.15lf x%d ", lp_obj_coef[c], c );
  fprintf ( file, ";\n" );
  
  for ( r = 0; r < lp_rows; r++ ) {
    for ( c = 0; c < lp_cols; c++ ) 
      if ( ! Equal( lp_coef_matrix[r][c], 0.0, LP_EPSILON ))
        fprintf( file, "%+.15lf x%d ", lp_coef_matrix[r][c], c );
    
    fprintf ( file, " %c %+.15lf;\n", lp_sense_out_char[r], lp_rhs[r] );

  } /* for r */

  fclose( file );

}  /* writeExternalLpSolve */
/**********************************************************************/
int solveExternalCplex(  ) {
  Abort( "solveExternalCplex() not yet implemented." );
} /* solveExternalCplex */
/**********************************************************************/
int solveExternalLpSolve(  ) {
/*
  Externally solves the LP in a file and then reads in the results.
  Returns either: LP_OPTIMAL, LP_INFEASIBLE or LP_UNBOUNDED.
*/
  FILE *file;
  char line[201];
  char cmd[80], str1[20], str2[20], str3[20], str4[20];
  int c, var;
  double val;
  int *state_seen;

  sprintf( cmd, "lp_solve/lp_solve < %s > %s", 
           LP_TEMP_FILENAME, SOLN_TEMP_FILENAME );
  system( cmd );

  if ((file = fopen( SOLN_TEMP_FILENAME, "r")) == NULL) 
    Abort( "Cannot open LP solution temp file for reading." );

  /* Whether solution or not, there are four strings, either:
     
     "Value of objective function:"

     or

     "This problem is infeasible"
     or

     "This problem is unbounded"

     One problem I see with the random generation of LPs with
     equalities is that you can get contradictory bounds. i.e.,

	C1 x2 > B1
	C2 x2 = B2

     Since the equality row has but one variable, it is completely
     determined.  When this happens, an empty solution file is
     generated, so we can catch these cases by saying that any empty
     solution file is infeasible.  We detect this by not being abole
     to read the 4 words.
     */

  if ( fscanf( file, "%s %s %s %s", str1, str2, str3, str4 ) != 4 ) {
    fclose( file );
    return( LP_INFEASIBLE );
  }

  if ( strcmp( str4, "infeasible" ) == 0 )  {
    fclose( file );
    return( LP_INFEASIBLE );
  }

  if ( strcmp( str4, "unbounded" ) == 0 )  {
    fclose( file );
    return( LP_UNBOUNDED );
  }

  fscanf( file, "%lf", &lp_obj_solution );

  state_seen = (int *) calloc( lp_cols, sizeof( int ) );

  while (fgets( line, 200, file) != NULL) {

    if ( sscanf( line, "x%d %lf\n", &var, &val ) <= 0 )
      continue;

    if ((var < 0)||(var >= lp_cols)) {
      fprintf(stderr, "Weird variable read from lp: %d\n", var);
      exit(-1);
    }
    if (state_seen[var] > 0 ) {
      fprintf(stderr, "Whoa, state %d seen twice!\n", var);
      exit(-1);
    }

    lp_solution[var] = val;
    state_seen[var]++;

  }

  fclose ( file );

  free ( state_seen );

  return ( LP_OPTIMAL );
}  /* solveExternalLpSolve */
/**********************************************************************/
void randomLPTest( int sparse ) {
  int i, c, external_status, internal_status;
  LP lp;
  int bad_lps = 0;
  int bad_result_match = 0;
  int bad_soln_match = 0;
  int num_infeasible = 0;
  int num_unbounded = 0;
  char msg[80];

  for ( i = 0; i < NUM_RANDOM_LPS; i++ ) {

    createRandomLP( sparse );

    if ( gAutoTestVerbose )
      printf( "[ Random LP size=%dx%d, sense=%d, has equal=%d ]\n",
              lp_rows, lp_cols, lp_obj_sense, lp_has_equality );
    
#ifdef HAVE_CPLEX
    writeExternalCplex();
#else
    writeExternalLpSolve();
#endif

#ifdef HAVE_CPLEX
    external_status = solveExternalCplex();
#else
    external_status = solveExternalLpSolve();
#endif

    lp = setFromRandomLP( sparse );
    
    internal_status = LP_solveLP( lp, NULL );

    /* Compare internal to external solution */

    if ( internal_status != external_status ) {

      bad_lps++;
      bad_result_match++;

      if ( gAutoTestVerbose )
        printf( "\tBad result match: E=%d, I=%d.\n",
                external_status, internal_status );

    } /* if LP return status does not match */

    else
      switch( internal_status ) {
      case LP_OPTIMAL:

        for ( c = 0; c < lp_cols; c++ )
          if ( ! Equal( lp->x[c], lp_solution[c], LP_SOLN_EPSILON )) {
            bad_lps++;
            bad_soln_match++;

            if ( gAutoTestVerbose ) {
              printf( "\tBad solution match.\n\tE: obj=%.2lf) ", 
                      lp->objval );
              for ( c = 0; c < lp_cols; c++ )
                printf( "%.2lf ", lp_solution[c] );
              printf( "\n" );
              
              printf( "\tI: obj=%.2lf) ", lp_obj_solution );
              for ( c = 0; c < lp_cols; c++ )
                printf( "%.2lf ", lp->x[c] );
              printf( "\n" );
            }
            
            break;
          }
        break;
      case LP_INFEASIBLE:
        num_infeasible++;
        break;
      case LP_UNBOUNDED:
        num_unbounded++;
        break;
      } /* switch */


    unlink( LP_TEMP_FILENAME );
    unlink( SOLN_TEMP_FILENAME );

    LP_freeLP( lp );

    freeRandomLP();

  } /* for i */
  
  if ( bad_lps > 0 ) {

    if ( gAutoTestVerbose ) {
      printf( "Num mismatched LP status = %d.\n", bad_result_match );
      printf( "Num mismatched LP solutions = %d.\n", bad_soln_match );
      printf( "Num infeasible LPs = %d.\n", num_infeasible );
      printf( "Num unbounded LPs = %d.\n", num_unbounded );
    }

    sprintf( msg,
             "There were %d LPs out of %d that gave different results.",
             bad_lps, NUM_RANDOM_LPS );
    error( msg );
  }

}  /* randomLPTest */
/**********************************************************************/
double *getRandomAlpha( double min, double max  ) {
   double *alpha;
   int j;

   alpha = newAlpha();

   for( j = 0; j < gNumStates; j++ )
     alpha[j] = getRandomDouble( min, max );

   return ( alpha );
}  /* *getRandomAlpha */
/**********************************************************************/
AlphaList getRandomAlphaList( int num_alphas, 
                              double min, double max  ) {
   AlphaList list = NULL;
   int i;
   
   list = newAlphaList();
   
   for( i = 0; i < num_alphas; i++ ) 
     prependAlphaList( list, getRandomAlpha(min,max), 0 );
   
   return( list );
}  /* getRandomAlphaList */
/**********************************************************************/
int encodeNumber( int z, int r ) {
/*
  When testing the neighbor and projection stuff, we will set the
  alpha vector component values so that by just looking at the value,
  we know which observation alpha list and which row the
  entry is in.  This allows us to make sure the computation is doing
  the proper thing.  We encode each row number with a specific set of
  bits and add one bit to allow us to differentiate between row zero
  and the absence of  soemthing from this observation.

  We divide an integer into a groups of bits.  Each group represents a
  particular observation.  The high order bit of the group is set just
  to indicate the presence of the observation (to differentiate the
  absence of a value for an observation from the presence of row 0 for
  an observation.)  The remaining bits are used to encode the row
  number. The total number of bits in a group is the number need to
  encode all the row numbers plus one. Examples (assume maximum of 8
  rows):

  Bit string             Obs     Row
  0000 0000 1000          0       0
  0000 0000 1011          0       3
  0000 0000 1100          0       4
  0000 1111 0000          1       7
  0000 1010 0000          1       2
  1101 0000 0000          2       5

  Thus, when a vector from a set of observations is selected and
  added, the component values of the vector can tell us if an
  observation vector set was used, and which row (vector) was used. 
  
*/
  int i;

  int num = 0x1;

  /* We first set the low order bits to reflect the presence of the
     observation and encode the row, then we shift the group of bits
     upward to the position they need to be in.

  /* This just sets the bit indicating the presence of an entry. */
  num = num << NEIGHBOR_ROW_BITS;

  /* This encode the row number. */
  num = num | r;

  /* This shifts the encoding over so it represents the particular
     observation. */
  num = num << ( z * (NEIGHBOR_ROW_BITS + 1));

  return ( num );
}  /* encodeNumber */
/**********************************************************************/
void decodeNumber( int num, int *present, int *row ) {
/*
  Uses the encoding scheme from encodeNumber().  The 'num' parameter
  is assumed to be the sum of encodings.

  Returns two arrays.  The present[z] array is booleans indicating if
  observation 'z' was used in the construction.  When present[z] is
  true, then row[z] indicates which row from that observations set was
  used. 
*/

  int obs_mask = 0x1;
  int row_mask, z;

  obs_mask = obs_mask << NEIGHBOR_ROW_BITS;
  row_mask = obs_mask - 1;

  for ( z = 0; z < gNumObservations; z++ ) {

    present[z] = num & obs_mask;
    row[z] = num & row_mask;

    num = num >> (NEIGHBOR_ROW_BITS + 1);

  } /* for z */

} /* decodeNumber */
/**********************************************************************/
AlphaList *makeEncodedProjection(  ) {
/*
  We use an encoding scheme for the first component of each vector,
  in each observation list, so we can track neighbors.  See
  encodeNumber() routine.
*/
  AlphaList *projection, node;
  int z, r;

  projection = (AlphaList *) calloc( gNumObservations,
                                     sizeof( *projection ));

  for ( z = 0; z < gNumObservations; z++ ) {
    
    projection[z] = getRandomAlphaList( NEIGHBOR_NUM_ROWS, 
                                        5.0, 10.0 );

    /* Encode the first component. */
    for ( r = 0, node = projection[z]->head; 
          node != NULL; r++, node = node->next )
      node->alpha[0] = encodeNumber( z, r );
      
  } /* for z */

  return ( projection );
}  /* *makeEncodedProjection */
/**********************************************************************/
AlphaList makeVector( int *row, AlphaList *projection ) {
/*
  Takes an array of row numbers (row[]) and constructs an alpha node
  that has that vector and the obs_source fields set.
*/
  AlphaList node, list;
  double *alpha;
  int i, z, r;

  alpha = (double *) calloc( gNumStates, sizeof(double));
  node = newAlphaNode( alpha, 0 );
  node->obs_source 
    = (AlphaList *) malloc( gNumObservations * sizeof( *node ));
  
  for ( z = 0; z < gNumObservations; z++ ) {

    r = 0;
    list = projection[z]->head;
    while (( r < row[z] ) && ( list != NULL )) {
      r++;
      list = list->next;
    } /* while */

    Assert( list != NULL, "Bad row in makeVector()." );
    
    node->obs_source[z] = list;

    for ( i = 0; i < gNumStates; i++ )
      alpha[i] += list->alpha[i];

  } /* for z */

  return ( node );

}  /* makeVector */
/**********************************************************************/
AlphaList makeEncodedAlphaList( int list_num ) {
  AlphaList list, node;
  double *alpha;
  int i;
  unsigned int n, value;
  
  list = newAlphaList();

  for ( n = 0; n < NUM_ENCODED_ALPHA_VECTORS; n++ ) {
    
    value = n <<  (NUM_ENCODED_SHIFT_BITS * list_num);

    alpha = newAlpha();
    for ( i = 0; i < gNumStates; i++ )
      alpha[i] = value + i/1000.0;

    node = newAlphaNodeObsSource( alpha, 0 );
    node->obs = list_num;
    appendNodeToAlphaList( list, node );

  } /* for n */

  return ( list );

}  /* makeEncodedAlphaList */
/**********************************************************************/
void matchPolicyGraph( AlphaList list, char *filename ) {
  FILE *file;
  int z, id, a, val;
  
  if ((file = fopen(filename , "r")) == NULL) {
    error( "The policy graph file: does not exist." );
    return;
  }
  
  list = list->head;
  while ( list != NULL ) {
    
    if ( list->obs_source == NULL ) {
      error( "NULL policy graph information found in node." );
      fclose( file );
      return;
    }
    
    /* If we run out of numbers before vectors, then this cannot be a
       match. */
    if ( fscanf( file, "%d %d", &id, &a ) < 2 ) {
      error( "Too many vectors." );
      fclose( file );
      return;
    }
    
    if (( id != list->id ) || ( a != list->action )) {
      error( "Non-matching ID or action." );
      fclose( file );
      return;
    }
    
    for ( z = 0; z < gNumObservations; z++ ) {

      /* A NULL pointer in the obs_source array means that the
         observation is not possible.  This is written out as an 'X'
         so we will skip over this one because we want to keep it
         simple and read only integers. */
      if ( list->obs_source[z] == NULL ) {

        while( fscanf( file, "%c", &val ) > 0 ) {

          /* Only expect to see spaces. Anything else is wrong or the
             'X' we are looking for. */
          if ( val != 32 )
            break;
        }

        if ( val != 'X' ) {
          error( "Didn't find 'X' for NULL choice." );
          fclose( file );
          return;
        }

        continue;
      } /* if obs_source[z] == NULL */

      if ( fscanf( file, "%d", &val ) < 1 ) {
        error( "Too many values in list." );
        fclose( file );
        return;
      }

      if ( val != list->obs_source[z]->id ) {
        error( "Non-matching observation choices." );
        fclose( file );
        return;
      }

    } /* for z */
    
    list = list->next;
  } /* while */

  /* Make sure we have exhausted the file too. */
  if ( fscanf( file, "%d", &val ) > 0 ) 
    error( "Not enough vectors found." );

  fclose( file );


}  /* matchPolicyGraph */
/**********************************************************************/

  /**********************************************************************/
 /**********    ALGORITHM TEST ROUTINES      ***************************/
/**********************************************************************/

/**********************************************************************/
void initPomdpSolveTest( PomdpSolveParams param ) {

  /* Open the POMDP file and initialize things. */
  initializePomdp( param->param_filename, 
                   param->impossible_obs_epsilon );

   switch ( param->method ) {
   case Enumeration:
     initEnumeration( );
     break;
   case Witness:
     initWitness();
     break;

#ifdef ZZZ
   case TwoPass:
     initTwoPass( );
     break;
   case LinearSupport:
     initLinSupport( );
     break;
   case IncrementalPruning:
     initIncPrune( );
     break;
#endif

   default:
     break;
   }  /* switch gMethod */

   initGlobal();
   initCommon();

}  /* initPomdpSolveTest */
/**********************************************************************/
void cleanUpPomdpSolveTest( PomdpSolveParams param ) {

  switch ( param->method ) {
  case Enumeration:
    cleanUpEnumeration( );
    break;
  case Witness:
    cleanUpWitness();
    break;
    
#ifdef ZZZ
  case TwoPass:
    cleanUpTwoPass( );
    break;
  case LinearSupport:
    cleanUpLinSupport( );
    break;
  case IncrementalPruning:
    cleanUpIncPrune( );
    break;
#endif
    
  default:
    break;
  }  /* switch gMethod */
  
  cleanUpCommon();
  cleanUpGlobal();
  cleanUpPomdp();
  
}  /* cleanUpPomdpSolveTest */
/**********************************************************************/
void setAlgorithmVariation( PomdpSolveParams param,
                            int variation, char *variation_name ) {
/*
  For each algorithm we want to try things with some different
  parameter settings.  This routine returns a parameter structure with
  the setting for the method and variation number sent in.
*/

  switch( param->method ) {

    /*************************/
  case Enumeration:
    strcpy( variation_name, "Enumeration" );
    break;

    /*************************/
  case TwoPass:
    strcpy( variation_name, "Two Pass" );
    break;

    /*************************/
  case LinearSupport:
    strcpy( variation_name, "Linear Support" );
    break;

    /*************************/
  case Witness:
    strcpy( variation_name, "Witness" );
    break;

    /*************************/
  case IncrementalPruning:
    switch( variation ) {
    case 0:
      strcpy( variation_name, "Normal Incremental Pruning" );
      param->ip_type = NormalIp;
      break;
    case 1:
      strcpy( variation_name, "Restricted Region IP" );
      param->ip_type = RestrictedRegionIp;
      break;
    case 2:
      strcpy( variation_name, "Generalized Incremental Pruning" );
      param->ip_type = GeneralizedIp;
      break;
    } /* switch variation */
    break;
    
    /*************************/
  default:
    Abort( "Unrecognized method." );
    break;
  } /* switch method */

}  /* setAlgorithmVariation */
/**********************************************************************/
void runPomdpSolveTest( PomdpSolveParams param ) {
/*
  Runs some number of POMDP solving iterations.

  Important parameters to be set in the param struct:

  param->param_filename - The POMDP file.
  param->initial_policy_filename - The initial policy to use of empty
                                   string for default.
  param->alpha_filename - A prefix for the solution files to compare
                          the computed answer against.
*/
  AlphaList next_alpha_list, real_alpha_list;
  int epoch;
  char soln_filename[MAX_FILENAME_LENGTH],
    msg[MAX_MSG_LENGTH];
  double start_time, stop_time;

  Assert( param->param_filename != NULL, 
          "No POMDP file specified." );

  initPomdpSolveTest( param );

   if ( param->initial_policy_filename[0] == NULL_CHAR ) 
     real_alpha_list = getDefaultInitialPolicy();
   else
     real_alpha_list 
       = readAlphaList( param->initial_policy_filename, -1 );

   Assert ( real_alpha_list != NULL,
            "Couldn't open initial policy file." );

   start_time = getSecs();

  /* We will test the first 20 epochs of the problem sent in. */
  for ( epoch = 1; epoch <= param->horizon; epoch++ ) {

    printf( "." );
    fflush( stdout );

    next_alpha_list = improveV( real_alpha_list, param );

    destroyAlphaList( real_alpha_list );

    sprintf( soln_filename, "%s%d", param->alpha_filename, epoch );
    real_alpha_list = readAlphaList( soln_filename, 0 );
    
    if ( ! similarAlphaList( real_alpha_list, 
                             next_alpha_list, 1e-9 )) {
      printf( "\n" );
      sprintf( msg, "Algorithm %s differs on %s at epoch %d.",
               method_str[param->method], 
               param->param_filename, epoch );
      error( msg );
      destroyAlphaList( next_alpha_list );
      break;
    } /* if found a bad match */

    destroyAlphaList( next_alpha_list );
    
  } /* for epoch */

  stop_time = getSecs();
  printf( "(%.2lf secs)\n", stop_time - start_time );

  cleanUpPomdpSolveTest( param );

} /* runPomdpSolveTest */
/**********************************************************************/
void runPolicyGraphTest( PomdpSolveParams param ) {
/*
  Makes sure the policy graph 'obs_source' fields get set properly
  during the solution process.
*/
  AlphaList next_alpha_list, real_alpha_list;
  char soln_filename[MAX_FILENAME_LENGTH],
    msg[MAX_MSG_LENGTH];

  Assert( param->param_filename != NULL, 
          "No POMDP file specified." );

  initPomdpSolveTest( param );
  
  Assert( param->initial_policy_filename[0] != NULL_CHAR,
          "No initial policy specified." ) ;
  
  real_alpha_list 
    = readAlphaList( param->initial_policy_filename, -1 );
  
  sortAlphaList( real_alpha_list );
  renumberAlphaList( real_alpha_list );
  
  Assert ( real_alpha_list != NULL,
           "Couldn't open initial policy file." );
  
  next_alpha_list = improveV( real_alpha_list, param );
  
  sortAlphaList( next_alpha_list );
  renumberAlphaList( real_alpha_list );
  
  if ( ! similarAlphaList( real_alpha_list, 
                           next_alpha_list, 1e-9 ))
    error( "Non-finitely transiency initial policy given?" );
  
  /* Check these against the precomputed versions. */
  matchPolicyGraph( next_alpha_list, param->pg_filename );

  destroyAlphaList( real_alpha_list );
  destroyAlphaList( next_alpha_list );
  cleanUpPomdpSolveTest( param );

}  /* runPolicyGraphTest */
/**********************************************************************/

  /**********************************************************************/
 /**********    MAIN TEST ROUTINES           ***************************/
/**********************************************************************/

/**********************************************************************/
void autoRandomTest(  ) {
  int *slots;
  int i, rand_int;
  double rand_double;

  /* Intialize the psuedo-random number generator. */
  randomize();
  
  /* Make a bunch of slots, put the random numbers in them, then check
     to make sure things are fairly evenly distributed. */
  slots = (int *) calloc( NUM_RANDOM_SLOTS, sizeof( int ));

    /******************************/
   /* Test of fran()             */
  /******************************/
  for ( i = 0; i < NUM_RANDOM_NUMBERS; i++ ) {

    rand_double = fran();

    if (( rand_double < 0.0 ) || ( rand_double > 1.0 )) {
      error( "fran() out of range number. Skipping test." );
      return;
    }

    /* Put it in a slot */
    slots[ (int) (NUM_RANDOM_SLOTS * rand_double) ]++;

  } /* for i */

  for ( i = 0; i < NUM_RANDOM_SLOTS; i++ ) {

    /* If less than 1/2 the expected value are there give an
       indication. */
    if ( slots[i] < NUM_RANDOM_NUMBERS/(NUM_RANDOM_SLOTS*2.0) ) {
      error( "Random numbers from fran() not looking so random." );
      break;
    } /* if badly distributes */

  } /* for i */
  
    /******************************/
   /* Test of getRandomDouble()  */
  /******************************/
  for ( i = 0; i < NUM_RANDOM_SLOTS; i++ ) 
    slots[i] = 0;

  for ( i = 0; i < NUM_RANDOM_NUMBERS; i++ ) {

    rand_double = getRandomDouble( RAND_DOUBLE_MIN, RAND_DOUBLE_MAX );
    
    if (( rand_double < RAND_DOUBLE_MIN ) 
        || ( rand_double > RAND_DOUBLE_MAX )) {
      error( "getRandomDouble() out of range number. " );
      break;
    } /* if out of range */

    rand_double = (rand_double - RAND_DOUBLE_MIN) 
      / ( RAND_DOUBLE_MAX - RAND_DOUBLE_MIN );

    Assert( ( rand_double >= 0.0 ) && ( rand_double < 1.0 ),
            "Conversion didn't do the right thing." );
      
    /* Put it in a slot */
    slots[ (int) (NUM_RANDOM_SLOTS * rand_double) ]++;

  } /* for i */

  for ( i = 0; i < NUM_RANDOM_SLOTS; i++ ) {

    /* If less than 1/2 the expected value are there give an
       indication. */
    if ( slots[i] < NUM_RANDOM_NUMBERS/(NUM_RANDOM_SLOTS*2.0) ) {
      error( "Random numbers from getRandomDouble() not looking so random." );
      break;
    } /* if badly distributes */

  } /* for i */
  
    /******************************/
   /* Test of getRandomInt()     */
  /******************************/
  for ( i = 0; i < NUM_RANDOM_SLOTS; i++ ) 
    slots[i] = 0;

  for ( i = 0; i < NUM_RANDOM_NUMBERS; i++ ) {

    rand_int = getRandomInt( RAND_INT_MIN, 
                             RAND_INT_MIN + NUM_RANDOM_SLOTS - 1 );
    
    if (( rand_int < RAND_INT_MIN ) 
        || ( rand_int > (RAND_INT_MIN + NUM_RANDOM_SLOTS - 1) )) {
      error( "getRandomInt() out of range number. " );
      break;
    } /* if out of range */

    /* Put it in a slot */
    slots[ rand_int - RAND_INT_MIN ]++;

  } /* for i */

  for ( i = 0; i < NUM_RANDOM_SLOTS; i++ ) {

    /* If less than 1/2 the expected value are there give an
       indication. */
    if ( slots[i] < NUM_RANDOM_NUMBERS/(NUM_RANDOM_SLOTS*2.0) ) {
      error( "Random numbers from getRandomInt() not looking so random." );
      break;
    } /* if badly distributes */

  } /* for i */

  free( slots );

}  /* autoRandomTest */
/**********************************************************************/
void autoPomdpTest(  ) {
  double val1;

  /****************************************/
  gValueType = REWARD_value_type;
  /****************************************/

  val1 = worstPossibleValue();
  if ( val1 > 0.0 )
    error( "worstPossibleValue() for rewards is wrong." );

  val1 = bestPossibleValue();
  if ( val1 < 0.0 )
    error( "bestPossibleValue() for rewards is wrong." );

  if ( isBetterValue( 1.0, 5.0, 1e-9 ))
    error( "isBetterValue() for rewards is wrong." );

  /****************************************/
  gValueType = COST_value_type;
  /****************************************/

  val1 = worstPossibleValue();
  if ( val1 < 0.0 )
    error( "worstPossibleValue() for costs is wrong." );

  val1 = bestPossibleValue();
  if ( val1 > 0.0 )
    error( "bestPossibleValue() for costs is wrong." );

  if ( isBetterValue( 5.0, 1.0, 1e-9 ))
    error( "isBetterValue() for costs is wrong." );

}  /* autoPomdpTest */
/**********************************************************************/
void autoAlphaTest(  ) {
/*
  Runs a bunch of test automaticaly checking that things work.
*/
  double *alpha1, *alpha2;
  int i;

  gNumStates = 10;
  gValueType = REWARD_value_type;

  alpha1 = newAlpha();
  alpha2 = newAlpha();
  for ( i = 0; i < gNumStates; i++ ) {
    alpha1[i] = i;
    alpha2[i] = i + 1;
  } 
  
  if ( ! sameAlpha( alpha1, alpha1, 1e-9 ))
    error( "sameAlpha() false on identical vectors." );

  if ( sameAlpha( alpha1, alpha2, 1e-9 ))
    error( "sameAlpha() true on different vectors." );
  
  free ( alpha1 );
  free ( alpha2 );

  alpha1 = newAlpha();
  for ( i = 0; i < gNumStates; i++ ) 
    alpha1[i] = i+1;
  alpha2 = duplicateAlpha( alpha1 );
  if ( ! sameAlpha( alpha1, alpha2, 1e-9 ))
    error( "duplicateAlpha() failed." );

  if ( isZeroAlpha( alpha1, 1e-9 ))
    error( "isZeroAlpha() true on non-zero vector." );
    
  for ( i = 0; i < gNumStates; i++ ) 
    alpha2[i] = 0.0;

  if ( ! isZeroAlpha( alpha2, 1e-9 ))
    error( "isZeroAlpha() false on zero vector." );
    
  if ( ! isDominatedVector( alpha1, alpha2 ))
    error( "isDominatedVector() false on better vector." );
  
  if ( isDominatedVector( alpha2, alpha1 ))
    error( "isDominatedVector() true on worse vector." );

  if ( ! isLexicographicallyBetterAlpha( alpha1, alpha2, 1e-9 ))
    error( "isLexicographicallyBetterAlpha() false on better vector." );
  
  if ( isLexicographicallyBetterAlpha( alpha2, alpha1, 1e-9 ))
    error( "isLexicographicallyBetterAlpha() true on worse vector." );

  free( alpha1 );
  free( alpha2 );
  
}  /* autoAlphaTest */
/**********************************************************************/
void autoAlphaListTest(  ) {
  AlphaList list1, list2, cur_node, node1, node2, node3;
  double *alpha1, *alpha2, *alpha3;
  double *b, val;
  int i, j;

  gNumStates = 10;
  gValueType = REWARD_value_type;

    /****************************************/
   /* Basic routines                       */
  /****************************************/

  alpha1 = newAlpha();
  alpha2 = newAlpha();
  alpha3 = newAlpha();

  for ( j = 0; j < gNumStates; j++ ) {
    alpha1[j] = 1.0;
    alpha2[j] = 2.0;
    alpha3[j] = 3.0;
  }

  node1 = newAlphaNode( alpha1, 0 );
  node2 = newAlphaNode( alpha2, 0 );
  node3 = newAlphaNode( alpha3, 0 );

  list1 = newAlphaList();

  appendNodeToAlphaList( list1, node1 );
  if ( list1->head != node1 || list1->tail != node1 )
    error( "appendNodeToAlphaList() on empty list failed." );

  appendNodeToAlphaList( list1, node2 );
  if ( list1->head == node2 || list1->tail != node2 )
    error( "appendNodeToAlphaList() on '1'-list failed." );

  prependNodeToAlphaList( list1, node3 );
  if ( list1->head != node3 || list1->tail == node3 )
    error( "preppendNodeToAlphaList() on list failed." );

  node1 = dequeueAlphaNode( list1 );
  if ( node1 != node3 )
    error( "dequeueAlphaNode() got wrong node." );

  enqueueAlphaNode( list1, node1 );
  if ( list1->tail != node1 )
    error( "enqueueAlphaNode() failed." );

  clearAlphaList( list1 );
  if ( list1->length != 0 || list1->head != NULL || list1->tail != NULL )
    error( "clearAlphaList() failed." );
    
  destroyAlphaList( list1 );

    /****************************************/
   /* Some of the more complex routines     */
  /****************************************/

  list1 = newIdentityAlphaList();
  b = (double *) calloc( gNumStates, sizeof( double ));
  b[4] = 1.0;
  val = bestVectorValuePrimed( list1, b, &node1, 
                               worstPossibleValue(), 1e-9 );
  if ( ! Equal( val, 1.0, 1e-9 ))
    error( "bestVectorValuePrimed() returned wrong value." );
  if ( !sameAlpha( node1->alpha, b, 1e-9 ))
    error( "bestVectorValuePrimed() returned wrong vector." );
  destroyAlphaList( list1 );
  free( b );

  list1 = newIdentityAlphaList();
  b = (double *) calloc( gNumStates, sizeof( double ));
  b[4] = 1.0;
  val = bestVectorValuePrimed( list1, b, &node1, 
                               bestPossibleValue(), 1e-9 );
  if ( val != bestPossibleValue() )
    error( "bestVectorValuePrimed() changed best value." );
  if ( node1 != NULL )
    error( "bestVectorValuePrimed() gives non-null on best value." );
  destroyAlphaList( list1 );
  free( b );
    
  list1 = newIdentityAlphaList();
  if ( ! sameAlphaList( list1, list1, 1e-9 ))
    error( "sameAlphaList() failed on same list." );

  list2 = getRandomAlphaList( gNumStates, 0.0, 1.0 );
  if ( ! similarAlphaList( list1, list1, 1e-9 ))
    error( "similarAlphaList() failed on same list." );
  if ( similarAlphaList( list1, list2, 1e-9 ))
    error( "similarAlphaList() succeeded on different lists." );
  destroyAlphaList( list2 );

  list2 = duplicateAlphaList( list1 );
  if ( ! sameAlphaList( list1, list2, 1e-9 ))
    error( "duplicateAlphaList() failed." );

  unionTwoAlphaLists( list1, list2 );
  if ( list1->length != (2 * gNumStates) )
    error( "unionTwoAlphaLists() failed." );
  destroyAlphaList( list1 );

  list1 = newIdentityAlphaList();
  alpha1 = newAlpha();
  for ( j = 0; j < gNumStates; j++ ) 
    alpha1[j] = 2.0;
  node1 = appendUniqueAlphaList( list1, alpha1, 0, 1e-9 );
  if ( node1 == NULL || list1->length != (gNumStates+1) )
    error( "appendUniqueAlphaList() failed for unique vector." );
  alpha2 = duplicateAlpha( alpha1 );
  node1 = appendUniqueAlphaList( list1, alpha2, 0, 1e-9 );
  if ( node1 != NULL || list1->length != (gNumStates+1)  )
    error( "appendUniqueAlphaList() failed for duplicate vector." );
  free( alpha2 );
  destroyAlphaList( list1 );
  
  list1 = newIdentityAlphaList();
  alpha1 = newAlpha();
  for ( j = 0; j < gNumStates; j++ ) 
    alpha1[j] = -9.0;
  if ( ! dominatedAlphaList( alpha1, list1 ))
    error( "dominatedAlphaList() failed for smaller vector." );
  for ( j = 0; j < gNumStates; j++ ) 
    alpha1[j] = 9.0;
  if ( dominatedAlphaList( alpha1, list1 ))
    error( "dominatedAlphaList() failed for larger vector." );
  free( alpha1 );
  destroyAlphaList( list1 );

  list1 = newIdentityAlphaList();
  clearMarkAlphaList( list1 );
  alpha1 = newAlpha();
  for ( j = 0; j < gNumStates; j++ ) 
    alpha1[j] = -9.0;
  if ( markDominatedAlphaList( alpha1, list1 ) != 0 )
    error( "markDominatedAlphaList() failed for smaller vector." );
  free( alpha1 );
  destroyAlphaList( list1 );

  list1 = newIdentityAlphaList();
  clearMarkAlphaList( list1 );
  alpha1 = newAlpha();
  for ( j = 0; j < gNumStates; j++ ) 
    alpha1[j] = 9.0;
  if ( markDominatedAlphaList( alpha1, list1 ) != gNumStates )
    error( "markDominatedAlphaList() failed for larger vector." );

  list2 = extractMarkedAlphaList( list1 );
  if ( list2->length != gNumStates || list1->length != 0 )
    error( "extractMarkedAlphaList() failed." );
  free ( alpha1 );
  destroyAlphaList ( list1 );
  destroyAlphaList ( list2 );

  list1 = newAlphaList();
  node1 = newAlphaNode( newAlpha(), 0 ); 
  cur_node = extractAlphaNode( list1, node1 );
  if ( cur_node != NULL )
    error( "extractAlphaNode() failed on empty list." );
  destroyAlphaNode( node1 );
  destroyAlphaList( list1 );

  list1 = newAlphaList();
  node1 = newAlphaNode( newAlpha(), 0 ); 
  appendNodeToAlphaList( list1, node1 );
  cur_node = extractAlphaNode( list1, node1 );
  if ( cur_node == NULL || list1->length != 0 )
    error( "extractAlphaNode() failed on single item list." );
  destroyAlphaNode( node1 );
  destroyAlphaList( list1 );

  list1 = createAlphaList( 3 );
  cur_node = extractAlphaNode( list1, list1->head );
  if ( cur_node == NULL || list1->length != 2 )
    error( "extractAlphaNode() failed for first item." );
  destroyAlphaNode( cur_node );
  destroyAlphaList( list1 );
  
  list1 = createAlphaList( 3 );
  cur_node = extractAlphaNode( list1, list1->head->next );
  if ( cur_node == NULL || list1->length != 2 )
    error( "extractAlphaNode() failed for middle item." );
  destroyAlphaNode( cur_node );
  destroyAlphaList( list1 );

  list1 = createAlphaList( 3 );
  cur_node = extractAlphaNode( list1, list1->tail );
  if ( cur_node == NULL || list1->length != 2 )
    error( "extractAlphaNode() failed for last item." );
  destroyAlphaNode( cur_node );
  destroyAlphaList( list1 );

  list1 = createAlphaList( 3 );
  node1 = newAlphaNode( newAlpha(), 0 );
  cur_node = extractAlphaNode( list1, node1 );
  if ( cur_node != NULL || list1->length != 3 )
    error( "extractAlphaNode() failed for non-existent node." );
  destroyAlphaNode( node1 );
  destroyAlphaList( list1 );


    /****************************************/
   /* Sorting                              */
  /****************************************/

  list1 = getRandomAlphaList( 25, 0.0, 10.0 );

  sortAlphaList( list1 );

  if ( list1->length != 25 )
    error( "sortAlphaList() lost or gained vectors." );
  
  list2 = list1->head;
  while ( list2 != NULL ) {
    
    if ( list2->next == NULL )
      break;

    if ( isLexicographicallyBetterAlpha( list2->alpha,
                                         list2->next->alpha,
                                         1e-9 )) {
      error( "sortAlphaList() didn't work." );
      break;
    }

    list2 = list2->next;
  } 

  destroyAlphaList( list1 );

} /* autoAlphaListTest */
/**********************************************************************/
void autoProjectionTest(  ) {
  AlphaList list;
  AlphaList **projection;
  AlphaList **file_list;
  int a, z;
  char filename[80];

  /* Don't need the POMDP until now. */
  initializePomdp( "examples/tiger.95.POMDP", 1e-9 );
  
  list = readAlphaList( "examples/tiger.95.optimal.alpha", -1 );

  if ( list == NULL ) {
    printf( "** Could not find initial alpha file. \n" );
    printf( "** Skipping test.\n" );
    return;
  }

  /* Call the make projection routine. */
  projection = makeAllProjections( list );

  /* Read saved projections from file. */
  file_list = allocateAllProjections();

  for ( a = 0; a < gNumActions; a++ )
    for ( z = 0; z < gNumObservations; z++ ) {
      
      sprintf( filename, "testing/tiger-a%d-z%d.alpha1", a, z );
      file_list[a][z] 
        = readAlphaList( filename, -1 );
      
      if ( file_list[a][z] == NULL ) {
        printf( "** Could not find projection[%d][%d] alpha file.\n",
                a, z );
        printf( "** Skipping test.\n" );
        return;
      }
    } /* for z */
  
  for ( a = 0; a < gNumActions; a++ )
    for ( z = 0; z < gNumObservations; z++ ) {

      if ( ! sameAlphaList( file_list[a][z], 
                            projection[a][z], 1e-9 ))
        error( "makeAllProjections() bad match to file." );
      
    } /* for z */
  
  freeAllProjections( file_list );
  freeAllProjections( projection );

  cleanUpPomdp();

}  /* autoProjectionTest */
/**********************************************************************/
void autoStatsTest(  ) {

  SolutionStats stat;

  stat = newSolutionStats( stdout, FALSE, FALSE );

  startContext( stat, 0 );
  recordLpStats( stat, 10, 4 );
  recordLpStats( stat, 10, 7 );
  recordLpStats( stat, 10, 13 );
  endContext( stat, 0 );

  startContext( stat, 1 );
  recordLpStats( stat, 10, 6 );
  recordLpStats( stat, 10, 3 );
  endContext( stat, 1 );

  startContext( stat, 2 );
  recordLpStats( stat, 10, 99 );
  endContext( stat, 2 );

  startContext( stat, 1 );
  recordLpStats( stat, 10, 87 );
  recordLpStats( stat, 10, 56 );
  recordLpStats( stat, 10, 1 );
  recordLpStats( stat, 10, 2 );
  recordLpStats( stat, 10, 3 );
  endContext( stat, 1 );

  startContext( stat, 0 );
  recordLpStats( stat, 10, 20 );
  endContext( stat, 0 );

  startContext( stat, 2 );
  recordLpStats( stat, 10, 2300 );
  recordLpStats( stat, 10, 2300 );
  recordLpStats( stat, 10, 45 );
  recordLpStats( stat, 10, 2100 );
  endContext( stat, 2 );

  if ( stat->lp_count[0] != 4
       || stat->constraint_count[0] != 44
       || stat->lp_count[1] != 7
       || stat->constraint_count[1] != 158
       || stat->lp_count[2] != 5
       || stat->constraint_count[2] != 6844 )
    error( "stats.c LP stat reporting problem." );

}  /* autoStatsTest */
/**********************************************************************/
void autoLpInterfaceTest(  ) {
  LP lp;
  int i;

    /***********************************/
   /* Dense Representation            */
  /***********************************/

  printf( "\tTesting dense LPs...\n" );
  randomLPTest( FALSE );

    /***********************************/
   /* Sparse Representation           */
  /***********************************/
  
  printf( "\tTesting sparse LPs...\n" );
  randomLPTest( TRUE );

}  /* autoLpInterfaceTest */
/**********************************************************************/
void autoRegionTest(  ) {
  AlphaList list;
  double *alpha;
  int i;
  PomdpSolveParams param;

  param = newPomdpSolveParams(  );

  gValueType = REWARD_value_type;
  gNumStates = REGION_TEST_STATES;

  /* We will create a set of N vectors of length N where the n'th
     vector has component 'n' set to 1.0 and all others set to
     zero. */
  list = newIdentityAlphaList();
  
  /* We create a few vectors that we know must have an empty
     region. */

  alpha = newAlpha();
  for ( i = 0; i < gNumStates; i++ )
    alpha[i] = 1.0 / ( 2.0 * gNumStates );

  if ( findRegionPoint( alpha, list, NULL, NULL, param ))
    error( "findRegionPoint() found false region point." );

  for ( i = 0; i < gNumStates; i++ )
    alpha[i] = 1.0 / ( gNumStates - 1.0 );

  if ( ! findRegionPoint( alpha, list, NULL, NULL, param ))
    error( "findRegionPoint() missed a region point." );

  free( alpha );
  destroyAlphaList( list );

}  /* autoRegionTest */
/**********************************************************************/
void autoParsimoniousTest(  ) {
  PomdpSolveParams param;
  AlphaList list, final_list, node, best_node;
  double *alpha, best_value;
  int i, j;

  gValueType = REWARD_value_type;
  gNumStates = PARSIMONIOUS_TEST_STATES;

    /***********************************/
   /*   isEmptyRegionSimpleCheck()    */
  /***********************************/

  /* Generate a set of random vectors. */
  list = getRandomAlphaList( PARSIMONIOUS_RAND_ALPHAS, 
                             -10.0, 15.0 );

  /* First select a random vector from the list and make sure
     simpleEmptyRegionCheck returns true. */
  for ( i = 0, node = list->head;
        i < getRandomInt( 0, PARSIMONIOUS_RAND_ALPHAS-1 ); 
        i++, node = node->next );

  if ( isEmptyRegionSimpleCheck( list, node->alpha, 1e-9, FALSE )
       == FALSE )
    error( "simpleEmptyRegionCheck() wrong result for vector in list." );

  /* Now make a vector that is better than all the ones in the list in
     one component to make sure the domination check doesn't return
     false positives. */
  alpha = newAlpha();
  for ( i = 0; i < gNumStates; i++ )
    alpha[i] = 0.0;

  /* Choose the component to use. */
  j = getRandomInt( 0, gNumStates-1 );
  alpha[j] = worstPossibleValue();
  node = list->head;
  while ( node != NULL ) {

    if ( isBetterValue( node->alpha[j], alpha[j], 1e-9 ))
      alpha[j] = node->alpha[j];
    
    node = node->next;
  } /* while */
  
  /* ensure it is larger by adding '1'. */
  alpha[j] += 1.0;

  if ( isEmptyRegionSimpleCheck( list, alpha, 1e-9, TRUE ) == TRUE )
    error( "simpleEmptyRegionCheck() true for non-empty region." );
  free( alpha );


  /* Now we make a vector that is component-wise worse than one of the
     vectors to make sure the domination checking takes care of this
     case. */

  alpha = newAlpha();
  /* First select a random vector from the list and make sure
     simpleEmptyRegionCheck returns true. */
  for ( i = 0, node = list->head;
        i < getRandomInt( 0, PARSIMONIOUS_RAND_ALPHAS-1 ); 
        i++, node = node->next );

  for ( i = 0; i < gNumStates; i++ )
    alpha[i] = node->alpha[i] - 1.0;

  if ( isEmptyRegionSimpleCheck( list, alpha, 1e-9, TRUE ) == FALSE )
    error( "simpleEmptyRegionCheck() false for empty region." );

  free( alpha );
  destroyAlphaList( list );

    /***********************************/
   /*   markBestAtSimplexVertices()   */
  /***********************************/

  /* Generate a set of vectors. */
  list = newIdentityAlphaList();
  clearMarkAlphaList( list );

  /* For the identity list, the simplex corners macth up dirrectly to
     the vectors.  We will save the "witness" point to iterate through
     them. */
  markBestAtSimplexVertices( list, TRUE, 1e-9 );
  
  for ( i = 0, node = list->head;
        node != NULL; 
        i++, node = node->next ) {
    
    /* All vectors should be marked. */
    if ( node->mark == FALSE ) {
      error( "Found missed vector in markBestAtSimplexVertices()" );
      break;
    }
    
    if ( node->witness == NULL ) {
      error( "No witness point set in markBestAtSimplexVertices()" );
      break;
    }
    
    for ( j = 0; j < gNumStates; j++ ) 
      if ( i == j ) {
        if ( ! Equal( node->witness[j], 1.0, 1e-9 ))
          error( "Wrong simplex set in markBestAtSimplexVertices()" );
      }
      else if ( ! Equal( node->witness[j], 0.0, 1e-9 ))
        error( "Wrong simplex set in markBestAtSimplexVertices()" );

    
  } /* for */ 

  destroyAlphaList( list );

    /***********************************/
   /*    markBestAtRandomPoints()     */
  /***********************************/

  /* Generate a set of random vectors. */
  list = newIdentityAlphaList();
  clearMarkAlphaList( list );

  /* First we use 1000 points to virtually ensure that each vector
     will be marked. */
  markBestAtRandomPoints( list, 1000, TRUE, 1e-9 );

  for ( node = list->head;
        node != NULL; 
        node = node->next ) {
    
    /* All vectors should most likely be marked (but not with
       probability 1.0). */
    if ( node->mark == FALSE ) {
      error( "Suspcious unmarked vector from markBestAtRandomPoints()" );
      break;
    }
    
    if ( node->witness == NULL ) {
      error( "No witness point set in markBestAtRandomPoints()" );
      break;
    }
    
    /* Make sure the witness point set really does yield this node as
       a vector. */
    best_node = bestVector( list, node->witness, &best_value, 1e-9 );
    
    if ( best_node != node )
      error( "Wrong witness point from markBestAtRandomPoints()" );
    
  } /* for */ 


  destroyAlphaList( list );

    /***********************************/
   /*       dominationCheck()         */
  /***********************************/

  /* Make a random alpha list and put a vector better than all of them
     in there and make sure we get only one vector back. */
  /* Generate a set of random vectors. */
  list = newIdentityAlphaList( );

  node = newAlphaNode( getRandomAlpha( 2.0, 5.0), 0 );

  appendNodeToAlphaList( list, node );

  dominationCheck( list ); 

  if (( list->length != 1 )
      || ( list->head != node )
      || (list->tail != node ))
    error( "DominationCheck() failed on list with one max." );

  destroyAlphaList( list );

  /* Now we put one in that is less than the rest and make sure we get
     one less vector. */

  list = newIdentityAlphaList();

  node = newAlphaNode( getRandomAlpha( -10.0, -5.0), 0 );

  appendNodeToAlphaList( list, node );

  dominationCheck( list ); 

  if ( list->length > gNumStates )
    error( "dominationCheck() failed to find one dom. vector." );

  if ( list->length < gNumStates )
    error( "dominationCheck() removed too many vectors." );

  /* Make sure the node is not in the list. */
  for ( node = list->head;
        node != NULL; 
        node = node->next ) 
    for ( i = 0; i < gNumStates; i++ )
      if ( node->alpha[i] < -1.0 )
       error( "dominationCheck() left bad vector in list." );
  
  destroyAlphaList( list );

    /***********************************/
   /*           prune()               */
  /***********************************/
  param = newPomdpSolveParams();
  param->use_witness_points = TRUE;
  param->prune_init_rand_points = 0;

  /* This first set of 50 vectors is parsimonious to begin with.  Make
     sure pruning does nothing to this set. */
  list = readAlphaList( "testing/prune-test1.alpha", 0 );
  Assert( list->length == 50, 
          "File 'prune-test1.alpha' supposed to have 50 useful vectors." );

  if ( prune( list, purge_prune, param ) != 0 )
    error( "prune() returns non-zero on parsimonious set." );

  if ( list->length != 50 )
    error( "prune() did not leave list unmodified." );

  destroyAlphaList( list );

  /* This next set of 100 vectors was randomly generated and there are
     31 useless vectors. Make sure that pruning leaves 91 vectors. */

  list = readAlphaList( "testing/prune-test-init.alpha", 0 );

  Assert( list->length == 100, 
          "Test alpha list not the expected size." );

  if ( prune( list, purge_prune, param ) != 9 )
    error( "prune() returns bad value on non-parsimonious set." );

  if ( list->length != 91 )
    error( "prune() left list with wrong number of vectors." );

  final_list = readAlphaList( "testing/prune-test-final.alpha", 0 );

  if ( ! similarAlphaList( list, final_list, 1e-9 ))
    error( "prune() left list not matching expected list." );

  destroyAlphaList( list );
  destroyAlphaList( final_list );

  /* This set of 100 vector has 99 useless ones.  Make sure we wind up
     with a single vector. */

  list = getRandomAlphaList( 99, 5.0, 10.0 );
  appendAlphaList( list, getRandomAlpha( 15.0, 20.0 ), 0 );

  if ( prune( list, purge_prune, param ) != 99 )
    error( "prune() did not remove the 99 useless vectors." );

  if ( list->length != 1 )
    error( "prune() left list with wrong number of vectors." );

  destroyAlphaList( list );
  free ( param );

}  /* autoParsimoniousTest */
/**********************************************************************/
void autoNeighborTest(  ) {
  PomdpSolveParams param;
  AlphaList list, node;
  AlphaList *projection;
  double *alpha;
  int i, z, r, *node_row, *neighbor_row, *neighbor_present, diff_count;

    /***********************************/
   /*   addNeighbor() routine         */
  /***********************************/

  gNumStates = 15;
  param = newPomdpSolveParams();
  param->domination_check = TRUE;

  list = getRandomAlphaList( 25, 5.0, 10.0 );

  /* Make sure adding a duplicate neighbor fails. */
  if ( addNeighbor( list, list->head->alpha,
                    NULL, 0, NULL, 
                    param->domination_check, 1e-9 ) == TRUE )
    error( "addNeighbor() allowed duplicate to be added." );

  /* Make sure adding a useless neighbor fails. */
  alpha = getRandomAlpha( 0.0, 4.5 );
  if ( addNeighbor( list, alpha,
                    NULL, 0, NULL, 
                    param->domination_check, 1e-9 ) == TRUE )
    error( "addNeighbor() allowed useless vector to be added." );
  free( alpha );

  /* Make sure adding a vector that is better than the rest results in
     an list of length 1. */
  alpha = getRandomAlpha( 20.0, 50.0 );
  if ( addNeighbor( list, alpha,
                    NULL, 0, NULL,
                    param->domination_check, 1e-9 ) == FALSE )
    error( "addNeighbor() denied adding a useful vector." );

  if ( sizeUnmarkedAlphaList( list ) != 1 )
    error( "addNeighbor() did not eliminate existing useless vectors." );

  destroyAlphaList( list );

  /* Make sure adding a vector works. The vector created is guaranteed
   to be useful at at least a point and useless at some other point. */
  list = getRandomAlphaList( 25, 5.0, 10.0 );
  alpha = newAlpha();
  alpha[0] = -50.0;
  alpha[1] = 50.0;
  if ( addNeighbor( list, alpha,
                    NULL, 0, NULL,
                    param->domination_check, 1e-9 ) == FALSE )
    error( "addNeighbor() denied adding a useful vector." );
  
  destroyAlphaList( list );

    /***********************************/
   /*   addAllNeighbor() routine      */
  /***********************************/

  /* We use an encoding scheme for the first component of each vector,
     in each observation list, so we can track neighbors. */
  gNumStates = 15;
  gNumObservations = 4;
  param->domination_check = FALSE;

  projection = makeEncodedProjection();

  list = newAlphaList();

  node_row = (int *) malloc( gNumObservations * sizeof( int ));

  for ( i = 0; i < gNumObservations; i++ )
    node_row[i] = getRandomInt( 0, NEIGHBOR_NUM_ROWS-1 );

  node = makeVector( node_row, projection );

  /* Because we are adding to an empty list, and we have
     gDominationCheck false, we know exactly how many vectors should
     be added. */
  if ( addAllNeighbors(list, node, projection, 
                       param->domination_check, 1e-9 ) 
       != gNumObservations * (NEIGHBOR_NUM_ROWS - 1))
    error( "Wrong number of neighbors added by addAllNeighbors()." );
  destroyAlphaNode( node );

  /* Now we make sure each neighbor generated is really a neighbor. */
  neighbor_row =  (int *) malloc( gNumObservations * sizeof( int ));
  neighbor_present =  (int *) malloc( gNumObservations * sizeof( int ));
  node = list->head;
  while ( node != NULL ) {

    /* A neighbor node will differ in only one position.  Use encoding
     to define this. */

    decodeNumber( node->alpha[0], neighbor_present, neighbor_row );

    diff_count = 0;
    for ( z = 0; z < gNumObservations; z++ ) {

      if ( ! neighbor_present[z] ) {
        error( "addAllNeighbors() - Missing vector for an observation." );
        break;
      }

      if ( neighbor_row[z] != node_row[z] )
        diff_count++;
      
    } /* for z */

    if ( diff_count != 1 )
      error( "addAllNeighbors() Bad neighbor added." );

    node = node->next;
  } /* while node != NULL */


  /* Clean up */
  destroyAlphaList( list );
  free( node_row );
  free( neighbor_row );
  free( neighbor_present );
  for ( z = 0; z < gNumObservations; z++ )
    destroyAlphaList( projection[z] );
  free( projection );

  free( param );
}  /* autoNeighborTest */
/**********************************************************************/
void autoCommonTest(  ) {
  AlphaList **projection, node, list;
  int a, z, m, i;
  double *alpha, *b, value;
  int num_vectors;
  
  extern int bestAlphaForBeliefQ( AlphaList node, double *b, 
                                  AlphaList *projection );
  
  /* Create a projection that has only one non-zero component per
     alpha vector.  */
  gNumObservations = 4;
  gNumActions = 4;
  num_vectors = 5;
  gNumStates = num_vectors + gNumActions;

  initCommon();

  projection = allocateAllProjections();

  /* Set up the projections in a special way so that we can set a
     belief state to easily determine what how the best vector should
     be constructed.
  */
  for ( a = 0; a < gNumActions; a++ ) 
    for ( z = 0; z < gNumObservations; z++ )  {
      projection[a][z] = newAlphaList();
      projection[a][z]->action = a;
      projection[a][z]->obs = z;

      for ( m = 0; m < num_vectors; m++ ) {
        alpha = newAlpha();
        
        for ( i = 0; i < gNumStates; i++ )
          if ( i == m )
            alpha[i] = 1.0 + i;
          else if ( i == (num_vectors + a ))
            alpha[i] = 1.0;
          else
            alpha[i] = 0.0;
        
        appendAlphaList( projection[a][z], alpha, a );

      } /* for m */
    } /* for z */
  
  b = (double *) calloc( gNumStates, sizeof( double ));
  node = newAlphaNode( newAlpha(), 0 );
  node->obs_source = newObsSourceArray();

    /***********************************/
   /*  bestAlphaForBeliefQ() routine  */
  /***********************************/

  for ( i = 0; i < gNumStates; i++ )
    b[i] = 1.0;

  for ( a = 0; a < gNumActions; a++ ) {

    if ( ! bestAlphaForBeliefQ( node, b, projection[a] ) )
      error( "bestAlphaForBeliefQ() bad return value." );
    
    for ( i = 0; i < gNumStates; i++ ) {
      
      if ( i == (num_vectors-1)) {
        
        if ( ! Equal( node->alpha[i], 
                      gNumObservations * num_vectors, 1e-9 ))
          error( "bestAlphaForBeliefQ() Alpha vector looks wrong." );
      }
      
      else if ( i == (num_vectors + a)) {

        if ( ! Equal( node->alpha[i], 
                      gNumObservations, 1e-9 ))
          error( "bestAlphaForBeliefQ() Alpha vector looks wrong." );
      }
      
      else if ( ! Equal( node->alpha[i], 0.0, 1e-9 )) {
        error( "bestAlphaForBeliefQ() Alpha vector looks wrong." );
        break;
      }

    } /* for i */

  } /* for a */

      /***********************************/
     /*     oneStepValue() routine      */
    /*             and                 */
   /*    makeAlphaVector() routine    */
  /***********************************/


  for ( i = 0; i < gNumStates; i++ )
    if ( i < num_vectors )
      b[i] = 1.0;
    else
      b[i] = 0.0;
  
  list = newAlphaList();

  for ( a = 0; a < gNumActions; a++ ) {
    b[num_vectors+a] = 1.0;

    value = oneStepValue( b, projection, &node, 1e-9 );
    if ( ! Equal( value,
                  (double) gNumObservations * num_vectors 
                  + gNumActions, 
                  1e-9 ))
      error( "oneStepValue() computed wrong best value." );
    
    if ( node->action != a )
      error( "oneStepValue() computed wrong best action." );

     if ( makeAlphaVector( list, projection, b, 1e-9 ) == NULL )
      error( "makeAlphaVector() failed when it shouldn't have." );
    
    b[num_vectors+a] = 0.0;
  } /* for a */

  /* Repeating this should  not add anything to list. */
  for ( a = 0; a < gNumActions; a++ ) {
    b[num_vectors+a] = 1.0;

    if ( makeAlphaVector( list, projection, b, 1e-9 ) != NULL )
      error( "makeAlphaVector() succeeded when it shouldn't have." );
    
    if ( list->length != gNumActions )
      error( "makeAlphaVector() left list the wrong size." );
    
    b[num_vectors+a] = 0.0;
  } /* for a */

  
  destroyAlphaList( list );

    /***********************************/
   /*     addVectorAtBeliefQ()        */
  /***********************************/

  /* Essentially just setBestAlphaForBeliefQ() so we'll do no testing
     for now. */ 

    /***********************************/
   /*     initWithSimplexCornersQ()   */
  /***********************************/
  
  list = newAlphaList();
  
  if ( initWithSimplexCornersQ( list, projection[0], TRUE, 1e-9 )
       != num_vectors )
    error( "initWithSimplexCornersQ() added wrong number of vectors." );

  destroyAlphaList( list );

    /***********************************/
   /* initWithRandomBeliefPointsQ()   */
  /***********************************/

  /* Pretty much the same as initWithSimplexCornersQ() so no explicit
     test cases yet. */

  freeAllProjections( projection );
  free( b );
  destroyAlphaNode( node );

  cleanUpCommon();

}  /* autoCommonTest */
/**********************************************************************/
void autoCrossSumTest(  ) {
  AlphaList list1, list2, result1, walk_ptr;
  int i, check[2][NUM_ENCODED_ALPHA_VECTORS];
  int list1_row, list2_row;

  gNumObservations = 4;
  gNumStates = 4;

  /* These lists are encoded so that from the result set we can
     determine which rows from each list the vector came from. */
  list1 = makeEncodedAlphaList( 0 );
  list2 = makeEncodedAlphaList( 1 );
  list1->action = 0;
  list2->action = 0;
  renumberAlphaList( list1 );
  renumberAlphaList( list2 );

  result1 = crossSum( list1, list2, FALSE );

  if ( result1->length != ( list1->length * list2->length )) {
    error( "crossSum() resulted in wrong number of vectors." );
    return;
  }

  /* Make sure each vector appears in the list the proper number of
     times by checking that the source pointers are set properly. */

  /* Clear out the counts for each vector in each list. */
  for ( i = 0; i < NUM_ENCODED_ALPHA_VECTORS; i++ )
    check[0][i] = check[1][i] = 0;

  walk_ptr = result1->head;
  while( walk_ptr != NULL ) {

    if (( walk_ptr->first_source == NULL )
        || ( walk_ptr->second_source == NULL ) ) {
      error("crossSum() NULL source pointer(s).");
      break;
    }
    
    check[0][walk_ptr->first_source->id]++;
    check[1][walk_ptr->second_source->id]++;

    walk_ptr = walk_ptr->next;
  } /* while walk_ptr != NULL */

  for ( i = 0; i < NUM_ENCODED_ALPHA_VECTORS; i++ )
    if (( check[0][i] != NUM_ENCODED_ALPHA_VECTORS)
        || ( check[1][i] != NUM_ENCODED_ALPHA_VECTORS ))
      error( "crossSum() source pointers not consistent." );

  /* For each vector we will compute the row from each list it seems
     to come from based on the values, then we match this to the
     source pointers. */

  /* Clear out the counts for each vector in each list. */
  for ( i = 0; i < NUM_ENCODED_ALPHA_VECTORS; i++ )
    check[0][i] = check[1][i] = 0;

  walk_ptr = result1->head;
  while( walk_ptr != NULL ) {
    
    list1_row
      = ((int) walk_ptr->alpha[0]) % NUM_ENCODED_ALPHA_VECTORS;
    list2_row
      = ((int) walk_ptr->alpha[0]) / NUM_ENCODED_ALPHA_VECTORS;

    if (( walk_ptr->first_source->id != list1_row )
        || ( walk_ptr->second_source->id != list2_row ))
      error( "crossSum() values don't match source settings." );

    check[0][list1_row]++;
    check[1][list2_row]++;

    walk_ptr = walk_ptr->next;
  } /* while walk_ptr != NULL */


  destroyAlphaList ( result1 );
  destroyAlphaList( list1 );
  destroyAlphaList( list2 );

}  /* autoCrossSumTest */
/**********************************************************************/
void autoPomdpSolveTest(  ) {
  PomdpSolveParams param;
  int i, v;
  char variation_name[80];

  param = newPomdpSolveParams();
  
  for ( i = 0; i < NUM_SOLVING_TESTS; i++ ) { 

    strcpy( param->param_filename, SolveParamStr[i] );
    strcpy( param->initial_policy_filename, SolveInitValStr[i] );
    strcpy( param->alpha_filename, SolveSolnStr[i] );

    printf( "\t%s\n", SolveDescriptStr[i] );

    for ( param->method = 0; 
          param->method < MAX_NUM_METHODS; 
          (param->method)++ ) {

      if ( ! AlgorithmTesting[param->method] ) {
        printf( "\t\t%s (skipping)\n", method_str[param->method] );
        continue;
      }
      
      param->horizon = SolveHorizon[i][param->method];

      for ( v = 0; v < NumAlgorithmVariations[param->method]; v++ ) {

        setAlgorithmVariation( param, v, variation_name );
        
        printf( "\t\t%s ", variation_name );
        
        runPomdpSolveTest( param );
      
      } /* for v */
      
    } /* for method */

  } /* for i */
  
}  /* autoPomdpSolveTest */
/**********************************************************************/
void autoPolicyGraphTest(  ) {

  PomdpSolveParams param;
  int i;

  param = newPomdpSolveParams();
  
  for ( i = 0; i < NUM_PG_TESTS; i++ ) { 

    strcpy( param->param_filename, PgParamStr[i] );
    strcpy( param->initial_policy_filename, PgInitValStr[i] );
    strcpy( param->pg_filename, PgSolnStr[i] );

    printf( "\t%s\n", PgDescriptStr[i] );

    for ( param->method = 0; 
          param->method < MAX_NUM_METHODS; 
          (param->method)++ ) {
      
      if ( ! AlgorithmTesting[param->method] ) {
        printf( "\t\t%s (skipping)\n", method_str[param->method] );
        continue;
      }
      
      printf( "\t\t%s\n", method_str[param->method] );

      runPolicyGraphTest( param );

    } /* for method */

  } /* for i */
  
}  /* autoPolicyGraphTest */
/**********************************************************************/
void initAutoTest(  ) {

  /* Throughout the course of testing we will alter the gNumStates
     variable.  Since we will use the temporary memory allocated in
     global.c we better make sure that we allocated enough space for
     testing. We just set it to the upper bound. */
  gNumStates = 200;
  gNumActions = 25;
  gNumObservations = 25;

  initGlobal();

}  /* initAutoTest */
/**********************************************************************/
void cleanUpAutoTest(  ) {

  gNumStates = 200;
  gNumActions = 25;
  gNumObservations = 25;

  cleanUpGlobal();

}  /* cleanUpAutoTest */
/**********************************************************************/
int main( int argc, char **argv ) {
  double start_user, start_system, stop_user, stop_system;
  int prev_errors;

  getSecsDetail( &start_user, &start_system );

  printf(  "       _________________\n");
  printf(  "     <<<   AUTO-TEST   >>>\n");
  printf(  "  ---------------------------\n");

  initAutoTest();

#ifndef UMMAGUMMA  

  printf( "Testing: random.c...\n" );
  prev_errors = gTotErrors;
  autoRandomTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );

  printf( "Testing: pomdp.c...\n" );
  prev_errors = gTotErrors;
  autoPomdpTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );

  printf( "Testing: alpha.c...\n" );
  prev_errors = gTotErrors;
  autoAlphaTest();
  autoAlphaListTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );

  printf( "Testing: projection.c...\n" );
  prev_errors = gTotErrors;
  autoProjectionTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );
  
  printf( "Testing: stats.c...\n" );
  prev_errors = gTotErrors;
  autoStatsTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );


  printf( "Testing: lp-interface.c...\n" );
  prev_errors = gTotErrors;
  autoLpInterfaceTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );

  printf( "Testing: region.c...\n" );
  prev_errors = gTotErrors;
  autoRegionTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );

  printf( "Testing: parsimonious.c...\n" );
  prev_errors = gTotErrors;
  autoParsimoniousTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );
  
  printf( "Testing: neighbor.c...\n" );
  prev_errors = gTotErrors;
  autoNeighborTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );
#endif
  
  printf( "Testing: common.c...\n" );
  prev_errors = gTotErrors;
  autoCommonTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );
  
  printf( "Testing: cross-sum.c...\n" );
  prev_errors = gTotErrors;
  autoCrossSumTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );
  
  /* The algorithm test will allocate their own memory, so clean up
     that which was used for the stuff above.  */
  cleanUpAutoTest(  );

  printf( "Testing: POMDP algorithms...\n" );
  prev_errors = gTotErrors;
  autoPomdpSolveTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );

  printf( "Testing: policy graph construction...\n" );
  prev_errors = gTotErrors;
  autoPolicyGraphTest();
  if ( gTotErrors == prev_errors )
    printf ( "\tNo errors.\n" );


  
  getSecsDetail( &stop_user, &stop_system );
  
  printf( "++++++++++++++++++++++++++++++++++++++++\n");
  printf( "Testing Finished. Total Errors = %d\n", gTotErrors );
  printf( "++++++++++++++++++++++++++++++++++++++++\n");

  reportTimes( stdout, stop_user - start_user,
              "User time =" );
  reportTimes( stdout, stop_system - start_system,
              "System time =" );
  reportTimes( stdout, 
               stop_user - start_user + stop_system - start_system,
                "Total execution time =" );

} /* main */
/**********************************************************************/

