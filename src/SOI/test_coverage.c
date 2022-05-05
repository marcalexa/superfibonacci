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
    FILE     *rotfile, *randfile;
    int     i,j;
    int        mode, numM, numRM=0;
    double   minminDist, maxminDist;
    double   *minDist, **M, **RM, in[5];
    char     inStr[120];
        
    if (argc < 2)
    {
        fprintf(stderr,"Usage:  test_coverage  rotfile  randfile  >  outfile \n");
        exit(1);
    }

    /* allocate data storage */
    
    numM = 0;
    rotfile = fopen(argv[1],"r");
    while (fgets(inStr,120,rotfile))
    {
    	numM++;
    }
    fclose(rotfile);
    
    mode = sscanf(inStr,"%lf %lf %lf %lf %lf",&in[0],&in[1],&in[2],&in[3],&in[4]);
    
    minDist = (double *) malloc(numM * sizeof(double));
    
    M = (double **) malloc(numM * sizeof(double *));
    for (i=0; i<numM; i++)
    {       
        M[i] = (double *) malloc(9 * sizeof(double));
    }
        
    numRM = 0;
    randfile = fopen(argv[2],"r");
    while (fgets(inStr,120, randfile))
    {
    	numRM++;
    }
    fclose(randfile);
    
    RM = (double **) malloc(numRM * sizeof(double *));   
    for (i=0; i<numRM; i++)
    {        
        RM[i] = (double *) malloc(9 * sizeof(double));
    }
       
    /* Read in the random rotations */
    
    numRM = LOAD_EA(argv[2],RM);
    
    /* Read in the sampled rotations */
    
    if (mode == 3)
    {        
        numM = LOAD_EA(argv[1],M);       
    }   
    else if (mode == 4)
    {       
        numM = LOAD_QUAT(argv[1],M);        
    }
    else  
    {
        fprintf(stderr,"Format problem:  use Euler angles or quaternions.  \n");
    }
        
    fprintf(stderr,"numRots: %i   numRand: %i\n",numM,numRM);
    
    COVERAGE(M, numM, RM, numRM, minDist);
    
    for (i=0;i<numRM;i++)
    {       
        printf("%15.6f \n",minDist[i]);
    }
          
    /* free memory */
    
    free(minDist);
    
    for (i=0; i<numM; i++)
    {        
        free(M[i]);
    }   
    
    free(M);
    
    for (i=0; i<numRM; i++)
    {       
        free(RM[i]);
    }
       
    free(RM);
    
    return 0;
}