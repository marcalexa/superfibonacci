/*****************************************************************************
* Author: Julie C. Mitchell
* Copyright (c) 2007 The University of Wisconsin
* All Rights Reserved
*
* Permission to use, copy, modify and distribute any part of this software for
* educational, research and non - profit purposes, without fee, and without a
* written agreement is hereby granted, provided that the above copyright
* notice, this paragraph and the following three paragraphs appear in all
* copies. In addition, the following reference must be cited in any pubications
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
*         Wisconsin Alumni Research Foundation (WARF)
*         614 Walnut Street, 13th floor
*         Madison, WI 53726
*         Telephone: 608 - 263 - 2500
*         Email: info@warf.org
*
* IN NO EVENT SHALL THE UNIVERSITY OF WISCONSIN BE LIABLE TO ANY PARTY FOR
* DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
* LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY
* OF WISCONSIN HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
* WISCONSIN HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
* ENHANCEMENTS, OR MODIFICATIONS. THE UNIVERSITY OF WISCONSIN MAKES NO
* REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR
* EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE
* SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
*
*****************************************************************************/

#include  <stdlib.h> 
#include  <stdio.h> 
#include  <math.h> 
#include "rotutils.h"

int main (int argc, char * argv[])
{
    /* variable defs */
    int i = 0;
    double  s, w[3];        /* current point from vector list */
    double  alpha;          /* increment angle for rotation sample */
    double  angles[3];      /* ordered angles */
    int n;         		    /* number of angular divisions */
    FILE *inFile;


    if (argc < 3)
    {
        fprintf(stderr,"Usage:  sample_soi   points.pts   orthogonal_divisions \n");
        exit(1);
    }

	inFile = fopen(argv[1],"r");
	n = atoi(argv[2]);
	if (n  <  1)
	{
		n = 1;
	}
	alpha = 2.0 * M_PI/((double) n);
	
    while(fscanf(inFile, "%lf %lf %lf", &w[0], &w[1], &w[2]) == 3)
    {
         
		if (w[2] > 1.0)
		{
			w[2] = 1.0;
		}
		
		if (w[2] < -1.0)
		{
			w[2] = -1.0;
		}

        angles[1] = acos(w[2]);
 
		// printf("%lf %lf %lf %lf\n",w[0],w[1],w[2],angles[1]);
        if (angles[1] == 0.0)
        {
			angles[2] = 0.0;
        }
        else
        {
			s = -w[0]/sin(angles[1]);
			if (s > 1.0)
			{
				s = 1.0;
			}
			
			if (s < -1.0)
			{
				s = -1.0;
			}
			
            angles[2] = acos(s);

            if(w[1] < 0)
            {
                angles[2] = 2.0 * M_PI - angles[2];
            }
        }
        if (angles[2]  >  M_PI)
        {
            angles[2] += - 2.0 * M_PI;
        }
       
        for (i=0; i<n; i++ )
        {
            angles[0] = - M_PI + i * alpha;
            fprintf(stdout, "%14.10lf\t %14.10lf\t %14.10lf\n", angles[0], angles[1], angles[2]);
        }
    }

    return 0;
} 
