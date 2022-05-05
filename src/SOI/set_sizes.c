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
#include  <string.h> 
#include  <math.h> 
#include "rotutils.h"

int num_icosamp(int L)
{
    return (12 + L*30 + 20*(L-1)*L/2);
}

int num_unifsamp(int n)
{
    return (floor(((double) n*n) / M_PI));
}



int main (int argc, char * argv[])
{
    /* variable defs */
    int i = 0;
	int n1, n2;
	int totnum = 0; 
	int done = 0;
	int L = 0;
	int N;

    if (argc < 4)
    {
        fprintf(stderr,"Usage:  set_sizes  -eis  angdivs  outname  >  run_soi.sh  \n");
        exit(1);
    }

	/* are we sampling by angle or by number of rotations? */
	
	if (strchr(argv[1],'n') != NULL)
	{
		totnum = 1;
	}
	
	/* convert angle or number of rotations to parameters for various sampling schemes */
	
	if (strchr(argv[1],'e') != NULL)
	{
		fprintf(stderr,"Laura and Roummel, please fill this in \n");
	}
	else if (strchr(argv[1],'i') != NULL)
	{
		if (totnum == 0)
		{
			while (!done)
			{
				N = floor(5.675 * (double) (L+1));
				if (N >= atoi(argv[2]))
				{
					done = 1;
				}
				else
				{
					L++;
				}
			}

			n1 = N;
			n2 = L;  
		}
		else
		{
			while (!done)
			{
				N = ((int) floor(5.675 * (double) (L+1))) * num_icosamp(L);
				if (N >= atoi(argv[2]))
				{
					done = 1;
				}
				else
				{
					L++;
				}
			}
			
			n1 = (int) floor(5.675 * (double) (L+1));
			n2 = L;
		}
		
		printf("./sample_icosa %i > %s.xyz\n",n2,argv[3]);
		printf("./sample_soi %s.xyz %i > %s.eul\n",argv[3],n1,argv[3]);
		
	}
	else if (strchr(argv[1],'s') != NULL)
	{
		if (totnum == 0)
		{
			n1 = atoi(argv[2]);
			n2 = num_unifsamp(n1);  
		}
		else
		{
			while (!done)
			{
				N = L * num_unifsamp(L);
				if (N >= atoi(argv[2]))
				{
					done = 1;
				}
				else
				{
					L++;
				}
			}
			
			n1 = L;
			n2 = num_unifsamp(n1);  
		}
		
		printf("./sample_sphere %i > %s.xyz\n",n2,argv[3]);
		printf("./sample_soi %s.xyz %i > %s.eul\n",argv[3],n1,argv[3]);
	}
	else
	{
		fprintf(stderr,"You must specify one of the flags -e, -i or -s.\n");
		exit(1);
	}
} 


