/*****************************************************************************
* Author: Julie C. Mitchell
* Copyright (c) 2007 The University of Wisconsin
* All Rights Reserved
* 
* Permission to use, copy, modify and distribute any part of this software for
* educational, research and non-profit purposes, without fee, and without a
* written agreement is hereby granted, provided that the above copyright
* notice, this paragraph and the following three paragraphs appear in all
* copies.  In addition, the following reference must be cited in any pubications 
* describing software or applications in which the codes have been used.
*
*
*
*
*
* 
* Those desiring to incorporate this software into commercial products or use
* for commercial purposes should contact 
*
*		Wisconsin Alumni Research Foundation (WARF)
*		614 Walnut Street, 13th floor
*		Madison, WI 53726
*		Telephone: 608-263-2500
*		Email: info@warf.org
* 
* IN NO EVENT SHALL THE UNIVERSITY OF WISCONSIN BE LIABLE TO ANY PARTY FOR
* DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
* LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY
* OF WISCONSIN HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
* WISCONSIN HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
* ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF WISCONSIN MAKES NO
* REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR
* EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE
* SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
*
**************************************************************************** */

#define ZERO 0.000000000000001

/* typedefs */
typedef double vec[3];
typedef int edge[2];
typedef int sim[3];

/* input and output */
int LOAD_EA(char *filename, double **M);
int LOAD_QUAT(char *filename, double **M);
int LOAD_PTS(char *filename, double **M);
void PRINTVEC(vec v);

/* convert input to matrices */
void QUAT_MAT(double *M, double *q);
void EA_MAT(double *M, double *E);
void AA_MAT(vec ax, double ang, vec *R);
void VV_MAT(vec vec1, vec vec2, vec *R);

/* matrix distances and metrics */
double MAT_DIST(double *M, double* N);
double PROJ_DIST(double *M, double* N);

/* matrix/vector manipulation */
void CROSS(vec v, vec w, vec vxw);
void MMULT (vec *Left, vec *Right, vec *R);
void ROTATE (vec *R, vec v);
void NORMAL(vec v);

/* tools for analyzing rotation samples */
void COVERAGE(double **M, int numM, double **RM, int numRM, double *minDist);
void SEPARATION(double **M, int numM, double *minDist, double *minminDist, double *maxminDist); 
void CALC_ENG(int numM, double **M, int numP, double *P, double *E) ;
 
