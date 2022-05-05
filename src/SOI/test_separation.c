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

/****************************************************************************
* Author: Julie Mitchell
* Function: ldist (main)
* Last Revision: 3/21/06
*
*  Calculates the least distance between each rotation and others in a sample
*
*  Usage:
*
*************************************************************************** */

int main(int argc, char **argv)
{    
    FILE     *myfile;
    int        i, j, mode, numM;
    double   minminDist,maxminDist;
    double   *minDist, **M, in[5];
    char	 inStr[120];
        
    if (argc < 2)
    {
        fprintf(stderr,"Usage:  test_separation  rotfile  >  outfile \n");
        exit(1);
    }


    numM = 0;
    
    /* allocate data storage */
    
    myfile = fopen(argv[1],"r");
    while (fgets(inStr,120,myfile))
    {
    	numM++;
    }
    fclose(myfile);
    
    mode = 2; //sscanf(inStr,"%lf %lf %lf %lf %lf",&in[0],&in[1],&in[2],&in[3],&in[4]);
    
    minDist = (double *) malloc(numM * sizeof(double));
    
    M = (double **) malloc(numM * sizeof(double *));
    for (i=0; i<numM; i++)
    {        
        M[i] = (double *) malloc(9 * sizeof(double));
        for (j=0; j<9; j++)
        {
        	M[i][j] = 0.0;
		}
    }
        
    /* Read in the sampled rotations */
    
    if (mode == 2)
    {
        /* points format */
        numM = LOAD_PTS(argv[1],M);
    }    
    else if (mode == 3)
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
        fprintf(stderr,"Format problem:  input points, Euler angles or quaternions.  \n");
    }
        
    fprintf(stderr,"numRots: %i\n",numM);
    
    /* calculate separation */
    
    SEPARATION(M,numM,minDist,&minminDist,&maxminDist);
        
    fprintf(stderr,"Minimum distance = %lf         Maximum minimum distance = %lf\n",minminDist,maxminDist);
    
    /* free memory */
    
    free(minDist);
    
    for (i=0; i<numM; i++)
    {
        free(M[i]);
    } 
    free(M);
    
    return 0;
}