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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <memory.h>
#include "rotutils.h"

int main(int argc, char **argv)
{
    FILE     *rotfile, *outfile;
    int     i,j,k,numM,numP;
    double    Energy, dist, p;
    int        mode;
    double   **M, *E, *P, in[5];
    char     inStr[120];
    
    if  (argc != 2)
    {
        printf("Usage:  test_energy   rotfile  >  outfile\n");
        return 0;
    }

    rotfile = fopen(argv[1], "r");

    numM = 0;
    
    /* allocate data storage */
    
    rotfile = fopen(argv[1],"r");
    while (fgets(inStr,120,rotfile))
    {
    	numM++;
    }
    fclose(rotfile);
    
    mode = sscanf(inStr,"%lf %lf %lf %lf %lf",&in[0],&in[1],&in[2],&in[3],&in[4]);
        
    M = (double **) malloc(numM * sizeof(double *));
    for (i=0; i<numM; i++)
    {
    	M[i] = (double *) malloc(9 * sizeof(double));
    }
    
    numP = 10;
    
    E = (double *) malloc(numP * sizeof(double));
    P = (double *) malloc(numP * sizeof(double));
    
    /* Read in the sampled rotations */
    
    if (mode == 3)
    {
        /* Euler angle format */
        numM = LOAD_EA(argv[1],M);        
    }    
    else if (mode == 4)
    {
        /* quaternion format */
        numM = LOAD_QUAT(argv[1],M);        
    }
    else  
    {
        fprintf(stderr,"Format problem:  use Euler angles or quaternions.  \n");
    }
        
    fprintf(stderr,"numRots: %i\n",numM);
    
    /* calculate repulsive energy for of various orders */
    
    P[0] = 0.25;           P[1] = 0.33;           P[2] = 0.5;
    P[3] = 1.0;            P[4] = 2.0;            P[5] = 3.0;
    P[6] = 4.0;            P[7] = 6.0;            P[8] = 9.0;
    P[9] = 12.0;
    
    CALC_ENG(numM, M,numP, P, E);
    
    printf("%15s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n","p=1/4","p=1/3","p=1/2","  p=1","  p=2","  p=3","  p=4","  p=6","  p=9"," p=12");
    printf("%15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n",E[0],E[1],E[2],E[3],E[4],E[5],E[6],E[7],E[8],E[9]);
    
    /* free memory */
    
    for (i=0; i<numM; i++)
    {        
        free(M[i]);        
    }
    free(M);    
    free(E);
    free(P);
    
    return 0;
}