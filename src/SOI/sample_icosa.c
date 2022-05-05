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
#include <math.h>
#include "rotutils.h"

void uvec(vec *iverts, int *thisface, vec bcen, vec v);
void set_icosa (vec *iverts, edge *iedges, sim *iface );
void anginterp(vec vec1, vec vec2, double a, vec p);

/****************************************************************************
* Author: Julie Mitchell
* Function: rot_sample
*
* Create a rotation sample with icosahedral symmetry.
*
* Variables:
*
* L        sampling density level relative to ligand mol
*
****************************************************************************/
int main (int argc, char **argv)
{
    
    
    int i, j, k, ii, jj, kk;
    int L, T;
    
    vec  iverts[12], vert1, vert2, bcen;
    sim  ifaces[20];
    edge iedges[30];
    
    double w, a, step, twist;
    vec p, p1, p2, R[3];
    int flat = 0;
    
    if (argc == 2)
    {
        L = atoi(argv[1]);
    }
    else
    {
        fprintf(stderr,"Usage:   sample_icosa  samplevel\n");
        exit(1);
    }

    /* define icosahedron */
    set_icosa(iverts, iedges, ifaces);
    
    /* loop over all vertices */
    for (i=0; i<12; i++)
    {
        /* print out results */
        PRINTVEC(iverts[i]);
        printf("\n");
    }
    
    /* loop over all edges */
    for (i=0; i<30; i++)
    {
        /* calculate subdivisions */
        for (j=1; j<L+1; j++)
        {
           /* take a linear combination of the two endpoints */
           
           a =  (double) j/L;
           
           bcen[0] = a;
           bcen[1] = 1.0-a;
           bcen[2] = 0.0;
            
           if (!flat)
            {
           		anginterp(iverts[iedges[i][0]],iverts[iedges[i][1]],a,p);   /* note: this could be made more efficient */  
			}
			else
			{                           
            	for (k=0; k<3; k++)
            	{
            	     p[k] = a * iverts[iedges[i][0]][k] + (1.0-a) * iverts[iedges[i][1]][k];
            	}
            }
			
            //PRINTVEC(iverts[iedges[i][0]]);
            //printf("\n");           
            //PRINTVEC(iverts[iedges[i][1]]);
            //printf("\n");           
            PRINTVEC(p);
            printf("\n");           
            //printf("\n");           
        }
    }
    
    /* loop over all simplices */
    for (i=0;i<20;i++)
    {
	
/*		w = 0.0;
		for (j=0;j<3;j++)
		{
			w+= iverts[ifaces[i][1]][j] * iverts[ifaces[i][0]][j];
		}
		w = 2.0 * M_PI / acos(w);
		fprintf(stderr,"%lf\n",w);
*/
		
        /* step along simplex edge using edge barystep j/L */
        for (j=1;j<L;j++)
        {
            a = (double) j/L;
            
            if (!flat)
            {
	           	anginterp(iverts[ifaces[i][1]],iverts[ifaces[i][0]],a,p1);       
	           	anginterp(iverts[ifaces[i][2]],iverts[ifaces[i][0]],a,p2);       
	            for (k=1;k<L-j;k++)
	            {
	                w = (double) k / (L-j);
	       	    	anginterp(p1,p2,w,p);  /* note: this could be made more efficient */ 
	                PRINTVEC(p);
	                printf("\n");
				}  
			} 
            else
            {
	            vert1[0] = 0.0;
	            vert1[1] = a;
	            vert1[2] = 1.0 - a;
	            
	            vert2[0] = 1.0 - a;
	            vert2[1] = a;
	            vert2[2] = 0.0;
	            
	            /* step across simplex interior using edge barystep k/(L-j) */
	            for (k=1;k<L-j;k++)
	            {
	                w = (double) k / (L-j);
	                
	                /* define Barycentric coords relative to simplex */
	                for (ii=0; ii<3; ii++)               
	                {
	                    bcen[ii] = w * vert1[ii] + (1.0 - w) * vert2[ii];                   
	                }
	                                
	                for (jj=0; jj<3; jj++)               
	                {
	                    p[jj] = 0.0;
	                    for (kk=0;kk<3;kk++)
	                    {
	                        p[jj] += bcen[kk] * iverts[ifaces[i][kk]][jj];
	                    }
	                }
	                /* print out results */
	                PRINTVEC(p);
	                printf("\n");
	           }
            }                 
         }
    }
    
    return 0;
}

           

/****************************************************************************
* Author: Julie Mitchell
* Function: anginterp
* Last Revision: 6/25/06
*
* Computes the unit vector defined by a point within a simplex of
* the icosahedron
*
* Variables:
*
* ivert vertices of the icosahedron
* asimp a simplex of the icosahedron
* bcen        barycentric coordinates for a point within the simplex
* v unit vector (output)
*
*************************************************************************** */

void anginterp(vec vec1, vec vec2, double a, vec p)

{    
  int i, j;
  vec ax, R[3];
  double angle;
    	
	CROSS(vec1, vec2, ax);
  NORMAL(ax);
	
	angle = 0.0;
	for (i=0; i<3; i++)
	{
		angle += vec1[i] * vec2[i];
	p[i] = vec1[i];
	}
	angle = a * acos(angle);
	
	if (fabs(angle-M_PI) < ZERO)
	{
		fprintf(stderr,"Warning: your matrix is not uniquely defined when start and end points are antipodal.\n");
	}

	AA_MAT(ax,angle,R);

  ROTATE(R,p);
  NORMAL(p);
           
  return;
}

/****************************************************************************
* Author: Julie Mitchell
* Function: uvec
* Last Revision: 6/25/06
*
* Computes the unit vector defined by a point within a simplex of
* the icosahedron
*
* Variables:
*
* ivert vertices of the icosahedron
* asimp a simplex of the icosahedron
* bcen        barycentric coordinates for a point within the simplex
* v unit vector (output)
*
*************************************************************************** */

void uvec(vec *iverts, int *thisface, vec bcen, vec v)
{    
    int i, j;
    
    for (i=0;i<3;i++)
    {       
        v[i] = 0.0;
        for (j=0;j<3;j++)
        {
        	v[i] += bcen[j] * iverts[thisface[j]][i];   
		}     
    }
       
    return;
}


/****************************************************************************
* Author: Julie Mitchell
* Function: set_icosa
* Last Revision: 6/25/06
*
* Define the vertices and simplices of an icosahedron
*
* Variables:
*
* ivert    vertices of the icosahedron (output)
* iface    indices defining the simplices of an icosahedron (output)
*
*************************************************************************** */

void set_icosa( vec *iverts, edge *iedges, sim *ifaces )
{       
    iverts[0][0] = -0.72361; iverts[0][1] = 0.52573;   iverts[0][2] = 0.44721;
    iverts[1][0] = -0.72361; iverts[1][1] = -0.52573;  iverts[1][2] = 0.44721;
    iverts[2][0] = -0.27639; iverts[2][1] = 0.85065;   iverts[2][2] = -0.44721;
    iverts[3][0] = -0.27639; iverts[3][1] = -0.85065;  iverts[3][2] = -0.44721;
    iverts[4][0] = 0.00000;  iverts[4][1] = 0.00000;   iverts[4][2] = 1.00000;
    iverts[5][0] = 0.00000;  iverts[5][1] = 0.00000;   iverts[5][2] = -1.00000;
    iverts[6][0] = 0.27639;  iverts[6][1] = 0.85065;   iverts[6][2] = 0.44721;
    iverts[7][0] = 0.27639;  iverts[7][1] = -0.85065;  iverts[7][2] = 0.44721;
    iverts[8][0] = 0.72361;  iverts[8][1] = 0.52573;   iverts[8][2] = -0.44721;
    iverts[9][0] = 0.72361;  iverts[9][1] = -0.52573;  iverts[9][2] = -0.44721;
    iverts[10][0] = 0.89443; iverts[10][1] = 0.00000;  iverts[10][2] = 0.44721;
    iverts[11][0] = -0.89443;iverts[11][1] = 0.00000;  iverts[11][2] = -0.44721;
    
    iedges[0][0] = 0;  iedges[0][1] = 1;
    iedges[1][0] = 0;  iedges[1][1] = 2;
    iedges[2][0] = 0;  iedges[2][1] = 4;
    iedges[3][0] = 0;  iedges[3][1] = 6;
    iedges[4][0] = 0;  iedges[4][1] = 11;
    iedges[5][0] = 1;  iedges[5][1] = 3;
    iedges[6][0] = 1;  iedges[6][1] = 4;
    iedges[7][0] = 1;  iedges[7][1] = 7;
    iedges[8][0] = 1;  iedges[8][1] = 11;
    iedges[9][0] = 2;  iedges[9][1] = 5;
    iedges[10][0] = 2; iedges[10][1] = 6;
    iedges[11][0] = 2; iedges[11][1] = 8;
    iedges[12][0] = 2; iedges[12][1] = 11;
    iedges[13][0] = 3; iedges[13][1] = 5;
    iedges[14][0] = 3; iedges[14][1] = 7;
    iedges[15][0] = 3; iedges[15][1] = 9;
    iedges[16][0] = 3; iedges[16][1] = 11;
    iedges[17][0] = 4; iedges[17][1] = 6;
    iedges[18][0] = 4; iedges[18][1] = 7;
    iedges[19][0] = 4; iedges[19][1] = 10;
    iedges[20][0] = 5; iedges[20][1] = 8;
    iedges[21][0] = 5; iedges[21][1] = 9;
    iedges[22][0] = 5; iedges[22][1] = 11;
    iedges[23][0] = 6; iedges[23][1] = 8;
    iedges[24][0] = 6; iedges[24][1] = 10;
    iedges[25][0] = 7; iedges[25][1] = 9;
    iedges[26][0] = 7; iedges[26][1] = 10;
    iedges[27][0] = 8; iedges[27][1] = 9;
    iedges[28][0] = 8; iedges[28][1] = 10;
    iedges[29][0] = 9; iedges[29][1] = 10;
    
    ifaces[0][0] = 0;  ifaces[0][1] = 1;  ifaces[0][2] = 4;
    ifaces[1][0] = 0;  ifaces[1][1] = 1;  ifaces[1][2] = 11;
    ifaces[2][0] = 0;  ifaces[2][1] = 2;  ifaces[2][2] = 6;
    ifaces[3][0] = 0;  ifaces[3][1] = 2;  ifaces[3][2] = 11;
    ifaces[4][0] = 0;  ifaces[4][1] = 4;  ifaces[4][2] = 6;
    ifaces[5][0] = 1;  ifaces[5][1] = 3;  ifaces[5][2] = 7;
    ifaces[6][0] = 1;  ifaces[6][1] = 3;  ifaces[6][2] = 11;
    ifaces[7][0] = 1;  ifaces[7][1] = 4;  ifaces[7][2] = 7;
    ifaces[8][0] = 2;  ifaces[8][1] = 5;  ifaces[8][2] = 8;
    ifaces[9][0] = 2;  ifaces[9][1] = 5;  ifaces[9][2] = 11;
    ifaces[10][0] = 2; ifaces[10][1] = 6; ifaces[10][2] = 8;
    ifaces[11][0] = 3; ifaces[11][1] = 5; ifaces[11][2] = 9;
    ifaces[12][0] = 3; ifaces[12][1] = 5; ifaces[12][2] = 11;
    ifaces[13][0] = 3; ifaces[13][1] = 7; ifaces[13][2] = 9;
    ifaces[14][0] = 4; ifaces[14][1] = 6; ifaces[14][2] = 10;
    ifaces[15][0] = 4; ifaces[15][1] = 7; ifaces[15][2] = 10;
    ifaces[16][0] = 5; ifaces[16][1] = 8; ifaces[16][2] = 9;
    ifaces[17][0] = 6; ifaces[17][1] = 8; ifaces[17][2] = 10;
    ifaces[18][0] = 7; ifaces[18][1] = 9; ifaces[18][2] = 10;
    ifaces[19][0] = 8; ifaces[19][1] = 9; ifaces[19][2] = 10;
    
    return;
}
