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
#include "rotutils.h"

int main(int argc, char **argv)
{
    FILE     *ptsfile, *engfile;
    double   P[3];
    double   eng, escale;
	double   xscale = 1.0;
	double   emin, emax;
    int      numP;
    

    if (argc < 2)
    {
        printf("Usage:  pts2pdb   points.pts   >   points.pdb \n");
        return 0;
    }

    ptsfile = fopen(argv[1], "r");
		
	if (argc > 2) 
	{
		xscale = atof(argv[2]);
	}
    
	if (argc > 3) 
	{
		engfile = fopen(argv[3],"r");
		
		emin =  10000.0;
		emax = -10000.0;
		
		while (fscanf(engfile,"%lf",&eng))
		{
			if (eng > emax)
			{
				emax = eng;
			}
			
			if (eng < emin)
			{
				emin = eng;
			}
		}
		
		rewind(engfile);
	}

    numP = 0;
       
    while (fscanf(ptsfile,"%lf %lf %lf",&P[0], &P[1], &P[2]) == 3)
    {
		if (argc > 3) 
		{
			fscanf(engfile,"%lf",&eng);
			eng = 99.0 * ( (eng - emin) / (emax - emin) );
		} 
		else
		{
			eng = P[0] * 49 + 50;
		}

        printf("ATOM  %5d  FE  SOI %5d    %8.3lf%8.3lf%8.3lf  1.00 %5.2lf\n", numP+1, numP+1, xscale * P[0], xscale * P[1], xscale * P[2], eng);
        numP++;
    }
    
    fclose(ptsfile);

	if (argc > 3) 
	{
		fclose(engfile);
    }
	
    return 0;
}