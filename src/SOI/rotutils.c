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
*        Wisconsin Alumni Research Foundation (WARF)
*        614 Walnut Street, 13th floor
*        Madison, WI 53726
*        Telephone: 608-263-2500
*        Email: info@warf.org
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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include "rotutils.h"

/****************************************************************************
* Author: Julie Mitchell
* Function: QUAT_MAT
* Last Revision: 3/21/07
*
*  Converts a quaternion, q, into a rotation matrix, M
*
*  Variables: M, a rotation matrix (1x9); q, a unit quaternion (1x4)
*
*************************************************************************** */
void QUAT_MAT(double* M, double *q)
{
    M[0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
    M[1] = 2*(q[1]*q[2] + q[0]*q[3]);
    M[2] = 2*(q[1]*q[3] - q[0]*q[2]);
    M[3] = 2*(q[1]*q[2] - q[0]*q[3]);
    M[4] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
    M[5] = 2*(q[2]*q[3] + q[0]*q[1]);
    M[6] = 2*(q[1]*q[3] + q[0]*q[2]);
    M[7] = 2*(q[2]*q[3] - q[0]*q[1]);
    M[8] = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
}

/****************************************************************************
* Author: Julie Mitchell
* Function: EA_MAT
* Last Revision: 3/21/07
*
*  Converts an Euler angle triplet, E, into a rotation matrix, M
*
*  Variables: M, a rotation matrix (1x9); E, an Euler angle triplet (1x3)
*
*************************************************************************** */
void EA_MAT(double* M, double *E)
{
    double phi, theta, psi;
    double cpsi,spsi,ctheta,stheta,cphi,sphi;
    phi = E[0];
    theta = E[1];
    psi = E[2];
    cpsi = cos(psi);
    spsi = sin(psi);
    ctheta = cos(theta);
    stheta = sin(theta);
    cphi = cos(phi);
    sphi = sin(phi);
    M[0] = (double)(  cpsi*cphi  -  spsi*ctheta*sphi  );
    M[1] = (double)( -cpsi*sphi  -  spsi*ctheta*cphi );
    M[2] = (double)(  spsi*stheta );
    M[3] = (double)(  spsi*cphi  +  cpsi*ctheta*sphi );
    M[4] = (double)( -spsi*sphi  +  cpsi*ctheta*cphi );
    M[5] = (double)( -cpsi*stheta );
    M[6] = (double)(  stheta*sphi );
    M[7] = (double)(  stheta*cphi );
    M[8] = (double)(  ctheta);
    return;
}

/****************************************************************************
 * Author: Julie Mitchell
 * Function: AA_MAT
 * Last Revision: 3/25/07
 *
 * Converts the axis-angle representation of a rotation to a matrix
 *
 * Variables: 
 *
 * ax	axis of rotation
 * ang	angle of rotation
 * R		rotation matrix (output) 
 * 
 *************************************************************************** */
void AA_MAT(vec ax, double ang, vec *R){

	double a, b, c;
	double CA, SA;

	a = ax[0]; b=ax[1]; c=ax[2];
	CA = cos(ang); SA = sin(ang);

	R[0][0] = a * a * (1.0 - CA) + CA;
	R[0][1] = a * b * (1.0 - CA) - c * SA;
	R[0][2] = a * c * (1.0 - CA) + b * SA;

	R[1][0] = a * b * (1.0 - CA) + c * SA;
	R[1][1] = b * b * (1.0 - CA) + CA;
	R[1][2] = b * c * (1.0 - CA) - a * SA;

	R[2][0] = a * c * (1.0 - CA) - b * SA;
	R[2][1] = b * c * (1.0 - CA) + a * SA;
	R[2][2] = c * c * (1.0 - CA) + CA;

return;}

/****************************************************************************
 * Author: Julie Mitchell
 * Function: VV_MAT
 * Last Revision: 3/25/07
 *
 * Determines that rotation that maps vec1 to vec2 along their orthogonal vector
 *
 * Variables: 
 *
 * vec1	first vector
 * vec2	second vector
 * R		rotation matrix (output) 
 * 
 *************************************************************************** */
void VV_MAT(vec vec1, vec vec2, vec *R){

	int i;
	vec ax;
	double angle;
	
	CROSS(vec1, vec2, ax);
	
	angle = 0.0;
	for (i=0; i<3; i++)
	{
		angle += vec1[i] * vec2[i];
	}
	angle = acos(angle);
	
	if (fabs(angle-M_PI) < ZERO)
	{
		fprintf(stderr,"Warning: your matrix is not uniquely defined when start and end points are antipodal.\n");
	}

	AA_MAT(ax,angle,R);

return;}

/****************************************************************************
* Author: Julie Mitchell
* Function: MAT_DIST
* Last Revision: 3/21/07
*
*  Computes the spherical distance between matrices M and N
*
*  Variables: M, a rotation matrix (1x9); N, a rotation matrix (1x9)
*
*************************************************************************** */
double MAT_DIST(double* M, double* N)
{
    double dist, cdist, tr;
    int i;
    tr = 0.0;
    for (i=0;i<9;i++)
    tr += M[i] * N[i];
    cdist = (tr - 1.0) / 2.0 ;
    if (cdist >   1.0) cdist =  1.0;
    if (cdist <  -1.0) cdist = -1.0;
    dist = acos(cdist);
    return dist;
}

/****************************************************************************
* Author: Julie Mitchell
* Function: PROJ_DIST
* Last Revision: 3/21/07
*
*  Computes the Euclidean distance between matrices M and N
*
*  Variables: M, a rotation matrix (1x9); N, a rotation matrix (1x9)
*
*************************************************************************** */
double PROJ_DIST(double* M, double* N)
{
    double dist, sqdist, dM;
    int i;
    for (i=0;i<9;i++)
    {
        dM = M[i] - N[i];
        sqdist += dM*dM;
    }
    dist = sqrt(sqdist);
    return dist;
}

/****************************************************************************
* Author: Julie Mitchell
* Function: LOAD_EA
* Last Revision: 3/22/07
*
*  Reads in a list of Euler angles and converts to rotation matrices
*
*  Variables: M, a set of rotation matrices (Nx9);
*
*************************************************************************** */
int LOAD_EA(char *filename, double **M)
{
    FILE *rotfile;
    double EA[3];
    int numM=0;
    rotfile = fopen(filename, "r");
    while (fscanf(rotfile,"%lf %lf %lf",&EA[0], &EA[1], &EA[2]) == 3)
    {
        EA_MAT(M[numM],EA);
        numM++;
    }
    fclose(rotfile);
    
    return numM;
}

/****************************************************************************
* Author: Julie Mitchell
* Function: LOAD_PTS
* Last Revision: 3/22/07
*
*  Reads in a list of points
*
*  Variables: M, a set of points (Nx3);
*
*************************************************************************** */
int LOAD_PTS(char *filename, double **M)
{
    FILE *ptsfile;
    int numP=0;

    ptsfile = fopen(filename, "r");
    while (fscanf(ptsfile,"%lf %lf %lf",&M[numP][0], &M[numP][1], &M[numP][2]) == 3)
    {
        numP++;
    }
    fclose(ptsfile);
    
    return numP;
}

/****************************************************************************
* Author: Julie Mitchell
* Function: LOAD_QUAT
* Last Revision: 3/22/07
*
*  Reads in a list of quaternions
*
*  Variables:
*
*************************************************************************** */
int LOAD_QUAT(char *filename, double **M)
{
    FILE *rotfile;
    double Q[4];
    int dummy, numM=0;
    rotfile = fopen(filename, "r");
    while (fscanf(rotfile,"%lf %lf %lf %lf",&Q[0], &Q[1], &Q[2], &Q[3]) == 4)
    {
        QUAT_MAT(M[numM],Q);
        numM++;
    }
    fclose(rotfile);
    return numM;
}

/****************************************************************************
* Author: Julie Mitchell
* Function: COVERAGE
* Last Revision: 3/22/07
*
*
*************************************************************************** */
void COVERAGE(double **M, int numM, double **RM, int numRM, double *minDist)
{
    int i,j;
    double dist, maxminDist;
    for (i=0;i<numRM;i++)
    {
        minDist[i] = 100.0;
        for (j=0;j<numM;j++)
        {
            dist = MAT_DIST(RM[i],M[j]);
            if (dist < minDist[i])
            minDist[i] = dist;
        }
    }
    /* calculate maximum minimum distance */
    maxminDist = -100.0;
    for (i=0;i<numRM;i++)
    {
        if (minDist[i] > maxminDist)
        maxminDist = minDist[i];
    }
    return;
}

/****************************************************************************
* Author: Julie Mitchell
* Function: SEPARATION
* Last Revision: 3/22/07
*
*
*************************************************************************** */
void SEPARATION(double **M, int numM, double *minDist, double *minminDist, double *maxminDist)
{
    int i,j;
    double dist;
    for (i=0;i<numM;i++)     /* initialize */
    minDist[i] = 100.0;
    for (i=0;i<numM;i++)
    {
        for (j=i+1;j<numM;j++)
        {
            /* calculate the pairwise matrix distances */
            dist = PROJ_DIST(M[i],M[j]);
            if (dist == 0)
            {
                fprintf(stderr,"Warning: duplicates at %i and %i (d=%f)\n",i,j,dist);
            }
            
            /* if this distance is less than the minimum, update */
            if (dist < minDist[i])
            {
            		minDist[i] = dist;
			}
			
            if (dist < minDist[j])
            {
            		minDist[j] = dist;
			}
        }
    }
    /* calculate maximum minimum distance */
    *minminDist =  100.0;
    *maxminDist = -100.0;
    for (i=0;i<numM;i++)
    {
        if (minDist[i] > *maxminDist)
        {
        	*maxminDist = minDist[i];
		}
		
        if (minDist[i] < *minminDist)
        {
        	*minminDist = minDist[i];
		}
    }
    return;
}

/****************************************************************************
* Author: Julie Mitchell
* Function: CALC_ENG
* Last Revision: 3/22/07
*
*  Calculate self repulsion energy of various orders
*
*  Variables:
*
*************************************************************************** */
void CALC_ENG(int numM, double **M, int numP, double *P, double *E)
{
    int i,j,k;
    double dist;
    for (k=0;k<numP;k++)
    E[k] = 0.0;
    for (i=0;i<numM;i++)
    {
        for (j=i+1;j<numM;j++)
        {
            dist = MAT_DIST(M[i],M[j]);
            dist = 2.0 * sin(dist/2.0);
            if (dist == 0)
            dist = 0.00001;
            for (k=0; k<numP; k++)
            E[k] += pow(dist,-P[k]);
        }
    }
    return;
}

/****************************************************************************
 * Author: Julie Mitchell
 * Function: PRINTVEC
 * Last Revision: 6/25/06
 *
 * Print a 3-component vector
 *
 * Variables: 
 *
 * v vector 
 * 
 *************************************************************************** */
 
void PRINTVEC (vec v){
 
        fprintf(stdout, "%15.10f %15.10f %15.10f ", (float) v[0], (float) v[1], (float) v[2]);

return;}

/****************************************************************************
 * Author: Julie Mitchell
 * Function: ROTATE
 * Last Revision: 6/25/06
 *
 * This function rotates a point according to a rotation matrix
 *
 * Variables:
 *
* R rotation matrix
 * v point in 3D space
 * 
 *************************************************************************** */

void ROTATE (vec *R, vec v) {
	
	int i, j;
	vec w;

	for (i=0; i<3; i++) {
		w[i] = 0.0;
		for (j=0; j<3; j++) w[i] += R[i][j]*v[j];
 	}
 
	for (i=0; i<3; i++) v[i] = w[i];
 
return;}

/****************************************************************************
 * Author: Julie Mitchell
 * Function: CROSS
 * Last Revision: 6/25/06
 *
 * This function computes a vector CROSS product.
 *
 * Variables: 
 *
 * v		the first unit vector
 * w		the second unit vector
 * vxw	the vector product of v with w 
 * 
 *************************************************************************** */
 
void CROSS(vec v, vec w, vec vxw) {

	vxw[0] = v[1]*w[2] - v[2]*w[1];
	vxw[1] = v[2]*w[0] - v[0]*w[2];
	vxw[2] = v[0]*w[1] - v[1]*w[0];

return;}


/****************************************************************************
* Author: Julie Mitchell
* Function: MMULT
* Last Revision: 6/25/06
*
* This function rotates a point according to a rotation matrix
*
* Variables:
*
* S matrix
* T matrix
* R		matrix
* 
*************************************************************************** */

void MMULT (vec *Left, vec *Right, vec *R) {

	int i, j, k;
	
	for (i=0;i<3;i++) {
		for (j=0;j<3;j++) R[i][j] = 0.0;		
		R[i][i] = 1.0;
		ROTATE(Right, R[i]);
		ROTATE(Left, R[i]);
	}
	
return;} 

/****************************************************************************
 * Author: Julie Mitchell
 * Function: NORMAL
 * Last Revision: 3/25/07
 *
 * Normalize a vector, so that it has unit length
 *
 * Variables: 
 *
 * v		a vector in 3D space
 * 
 *************************************************************************** */
void NORMAL(vec v) { 
	
	double vnorm; 
	int i;
	 
	vnorm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

	if (vnorm > 0.0)
	{
		for (i=0; i<3; i++) 
		{
			v[i] /= vnorm;
		}
	}	
	
	return;
}
