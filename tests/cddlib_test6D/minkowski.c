/* minkowski.c: Main test program to call the cdd library cddlib
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.93, July 15, 2003
   Standard ftp site: ftp.ifor.math.ethz.ch, Directory: pub/fukuda/cdd
*/

/*  This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include "setoper.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

int main() {

	void callCdd();

	clock_t start = clock();
	int i;
	for (i=0; i<1;i++){
		callCdd();
	}
	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf ("seconds: %f  \n", seconds);
	return 0;
}

void callCdd()
{
	dd_MatrixPtr M2=NULL, A=NULL, G=NULL;
	dd_PolyhedraPtr poly;
	dd_LPPtr lp;
  dd_ErrorType err=dd_NoError;
  FILE *reading=NULL, *writing;
  mytype val;
  time_t starttime,endtime;
  dd_rowset redrows,linrows;

  dd_set_global_constants();  /* First, this must be called. */

  	  dd_NumberType numb = dd_Real;   /* set a number type */
  	  dd_rowrange m=5;    /* number of rows, number of edges  */
  	  dd_colrange n=6;    /* number of columns, dimension of the space */
  	  dd_MatrixPtr M=dd_CreateMatrix(m,n);

  	  dd_set_si(M->matrix[0][0], (long)1.0);
  	  dd_set_si(M->matrix[0][1], (long)1.0);
  	  dd_set_si(M->matrix[0][2], (long)1.0);
  	  dd_set_si(M->matrix[0][3], (long)0.0);
  	  dd_set_si(M->matrix[0][4], (long)0.0);
  	  dd_set_si(M->matrix[0][5], (long)0.0);

  	  dd_set_si(M->matrix[1][0], (long)-1.0);
  	  dd_set_si(M->matrix[1][1], (long)1.0);
  	  dd_set_si(M->matrix[1][2], (long)1.0);
  	  dd_set_si(M->matrix[1][3], (long)0.0);
  	  dd_set_si(M->matrix[1][4], (long)0.0);
  	  dd_set_si(M->matrix[1][5], (long)0.0);

  	  dd_set_si(M->matrix[2][0], (long)1.0);
  	  dd_set_si(M->matrix[2][1], (long)-1.0);
  	  dd_set_si(M->matrix[2][2], (long)1.0);
  	  dd_set_si(M->matrix[2][3], (long)0.0);
  	  dd_set_si(M->matrix[2][4], (long)0.0);
  	  dd_set_si(M->matrix[2][5], (long)0.0);

  	  dd_set_si(M->matrix[3][0], (long)-1.0);
  	  dd_set_si(M->matrix[3][1], (long)-1.0);
  	  dd_set_si(M->matrix[3][2], (long)1.0);
  	  dd_set_si(M->matrix[3][3], (long)0.0);
  	  dd_set_si(M->matrix[3][4], (long)0.0);
  	  dd_set_si(M->matrix[3][5], (long)0.0);

  	  dd_set_si(M->matrix[4][0], (long)0.0);
  	  dd_set_si(M->matrix[4][1], (long)0.0);
  	  dd_set_si(M->matrix[4][2], (long)1.0);
  	  dd_set_si(M->matrix[4][3], (long)0.0);
  	  dd_set_si(M->matrix[4][4], (long)0.0);
  	  dd_set_si(M->matrix[4][5], (long)0.0);

  	  M->representation=dd_Generator;
  	  poly=dd_DDMatrix2Poly(M, &err);  /* compute the second (generator) representation */
	  dd_init(val);

	  dd_rowindex newpos;
	  dd_rowset impl_lin,redset;

	  dd_MatrixCanonicalize(&M, &impl_lin, &redset, &newpos, &err);
	  dd_WriteMatrix(stdout, M);

//	  set_addelem(M->linset,1); /* setting the first to be equality */

//	  dd_DDInputAppend(&poly,M, &err); /* append the two inequalities and compute the generators */
	  if (err!=dd_NoError) goto _L99;
	  A=dd_CopyInequalities(poly);  /* get the inequalities (=input). */
	  G=dd_CopyGenerators(poly);  /* get the generators (=output). */
	  printf("\nNew H-representation with added inequalities:\n");
	  dd_WriteMatrix(stdout,A);  printf("\n");
	  dd_WriteMatrix(stdout,G);

	  _L99:
	    if (err!=dd_NoError) dd_WriteErrorMessages(stdout,err);
}


/* end of simplecdd.c */
