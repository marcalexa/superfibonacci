/*****************************************************************************
* Authors: Laura LeGault, Roummel F. Marcia and Julie C. Mitchell
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
*               Wisconsin Alumni Research Foundation (WARF)
*               614 Walnut Street, 13th floor
*               Madison, WI 53726
*               Telephone: 608-263-2500
*               Email: info@warf.org
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
#include <stdlib.h>
#include <math.h>
#include <memory.h>

double   n_choose_k( int n, int k);
double   integrate_sin_2n(int n, double theta0, double theta1);
double   arclength( double a, double b, double theta0, double theta1);
double   compute_theta( double a, double b, double theta_start, double arclength_start, double theta0, double length_targ);
void     get_angles_slice(double thetas[], int nrings, double a, double b, double h);
double   interp_quad( double t1, double t2, double t3, double r1, double r2, double r3, double al);
void     get_angles_ellipse( double a, double b, double thetas[], int stagger, double nodes);
void     order_abc(double x, double y, double z, double *a, double *b, double *c, int index[]);

int main(int argc, char *argv[])
{
    /* variables */
    double   a, b, c, r, p_sa;
    double   a1, b1, c1;
    int      aind[3];
    double   alpha = 2.0/3.0;

    int      n, n_tri, n_ring, total, num, npts;
    double   h, t, x, y, z;
    double   s1, s2, s3, p1, p2, p3;
    double   area, circ, cir, side, hr, denom, hslice, hr0, denom0, circ0;
    double   flat_ell, majax, minax;
    double   **pts;
    int      i, j, k;
    
    if(argc != 5)
    {
        printf("Usage:  sample_ellipsoid   x_length   y_length   z_length   number_of_points\n");
        return 0;
    }
    
    a1 = atof(argv[1]);
    b1 = atof(argv[2]);
    c1 = atof(argv[3]);
    n  = atoi(argv[4]);
    
    order_abc(a1, b1, c1, &a, &b, &c, aind);

    if (n < 5)
    {
        printf("Number of nodes must be greater than four.\n");
        exit(1);
    }
    
    if (a <= 0.0 || b <= 0.0 || c <= 0.0)
    {
        printf("Lengths of principal axes must be positive.\n");
        exit(1);
    }
    
    r  = 1.6075;
    s1 = pow(a,r);
    s2 = pow(b,r);
    s3 = pow(c,r);
    p1 = s1*s2;
    p2 = s1*s3;
    p3 = s2*s3;
    p_sa   = (p1+p2+p3)/3.0;
    p_sa   = pow(p_sa, (1.0/r));
    p_sa   = (4.0*M_PI)*p_sa;
    n_tri  = 2*n - 4;
    area   = p_sa / (double) n_tri;
    h      = sqrt(sqrt(3.0)*area);
    side   = 2.0 * h / sqrt(3.0);
    hr     = pow( a-c, 2.0) / pow (a+c, 2.0);           // note c is used and not b
    denom  = 10.0 + sqrt( 4.0-(3.0*hr) );
    circ   = M_PI * (a + c) * (1.0 + (3.0*hr/denom));   // Ramanujan II for an ellipse
    circ   = circ/2.0;
    n_ring = (int) round(circ/h);

    flat_ell = b - (alpha*a + (1.0 - alpha)*sqrt(a*c));

    if (flat_ell > 0.0){
           hr0    = pow(b-c, 2.0) / pow (b+c, 2.0);
	   denom0 = 10.0 + sqrt( 4.0-(3.0*hr0) );
	   circ0  = M_PI * (b + c) * (1.0 + (3.0*hr0/denom0))/2.0;
	   n_ring = round(circ0/h);
    }

    double rings[n_ring];
    double nodes[n_ring];
    double thetas[n_ring];
    total  = 2;
    num    = n_ring + 1;
    t      = M_PI/num;
    hslice = circ/(double) num;


    if (flat_ell > 0.0)     get_angles_slice(thetas, n_ring, c, a, hslice);
    else                    get_angles_slice(thetas, n_ring, a, c, hslice);
           
    for (i = 0; i < n_ring; i++) 
    {

         if (flat_ell > 0.0){
	        rings[i] = c*cos(thetas[i]);
		x        = a*sqrt((1.0-rings[i]/c)*(1.0+rings[i]/c));
		y        = b*sqrt((1.0-rings[i]/c)*(1.0+rings[i]/c));
		hr       = pow( x-y, 2.0) / pow (x+y, 2.0);
		denom    = 10.0 + sqrt( 4.0-(3.0*hr) );
		cir      = M_PI * (x + y) * (1.0 + (3.0*hr/denom));   // Ramanujan II for an ellipse
	 } else {
	        rings[i] = a*cos(thetas[i]);
		y        = b*sqrt((1.0-rings[i]/a)*(1.0+rings[i]/a));
		z        = c*sqrt((1.0-rings[i]/a)*(1.0+rings[i]/a));
		hr       = pow( y-z, 2.0) / pow (y+z, 2.0);
		denom    = 10.0 + sqrt( 4.0-(3.0*hr) );
		cir      = M_PI * (y + z) * (1.0 + (3.0*hr/denom));   // Ramanujan II for an ellipse
	 }
	 nodes[i] = round(cir/side);
	 total    = total + (int) nodes[i];

    }

    npts = total;
    total--;

    pts = calloc( npts, sizeof(double));
    for (j = 0; j < npts; j++)     pts[j] = calloc( 3, sizeof(double));

    int    t_index = 0;    // must be set to 0 initially
    double *theta_ell;
    for (j = 0; j < n_ring; j++) 
    {
           if (flat_ell > 0.0)
	   {
	          majax = sqrt(pow(a,2.0)*(1.0 - pow(rings[j],2.0)/pow(c,2.0)));
		  minax = sqrt(pow(b,2.0)*(1.0 - pow(rings[j],2.0)/pow(c,2.0)));
        
		  // stagger the nodes                 
		  if (j > 0)     t_index = (t_index + 1)%2;
		  else 	         if ( ((int) nodes[j])%4 == 0)  t_index = 1;

	   } else {
	          majax = sqrt(pow(b,2)*(1 - pow(rings[j],2)/pow(a,2)));
		  minax = sqrt(pow(c,2)*(1 - pow(rings[j],2)/pow(a,2)));
                        
		  // stagger the nodes
		  if (j > 0){  
		        if ( (nodes[j] != nodes[j-1]) && (( (int) nodes[j]- (int) nodes[j-1])%2 == 0))      t_index = t_index;
			else        	                                                                    t_index = (t_index + 1)%2;
		  } else {  
		        if ( ((int) nodes[j])%4 == 0)  t_index = 1;
		  }
	   }
	   theta_ell = calloc( (int) nodes[j], sizeof(double));
	   get_angles_ellipse( majax, minax, theta_ell, t_index, nodes[j]);
		  
	   for (k = 0; k < nodes[j]; k++) 
	   {
	           if (flat_ell > 0.0)
		   {
		        pts[total][aind[0]] = majax*cos(theta_ell[k]);                       
			pts[total][aind[1]] = minax*sin(theta_ell[k]);
			pts[total][aind[2]] = rings[j];
		   } else {
		        pts[total][aind[0]] = rings[j];
			pts[total][aind[1]] = majax*cos(theta_ell[k]);
			pts[total][aind[2]] = minax*sin(theta_ell[k]);                  
		   }
		   total--;
		   if (total < 0)
		   {
			 printf("total = %d is less than total nodes!\n", total);
			 return 0;
		   }
	   }
	   free (theta_ell);
    }

    if (flat_ell > 0.0){
	   pts[total][aind[0]] = 0.0;
	   pts[total][aind[1]] = 0.0;
	   pts[total][aind[2]] =   c;
    } else {
           pts[total][aind[0]] = a;
	   pts[total][aind[1]] = 0.0;
	   pts[total][aind[2]] = 0.0;

    }
    total--;
    if (total != 0){
         printf("total = %d is not equal to total nodes!\n", total);
	 return 0;
    }
	
    if (flat_ell > 0.0){
	   pts[total][aind[0]] = 0.0;
	   pts[total][aind[1]] = 0.0;
	   pts[total][aind[2]] =  -c;
    } else {
           pts[total][aind[0]] = -a;
	   pts[total][aind[1]] = 0.0;
	   pts[total][aind[2]] = 0.0;
    }


    for (i = 0; i < npts; i++)
    {
        printf("%15.10lf %15.10lf %15.10lf\n",pts[i][0],pts[i][1],pts[i][2]);
    }


    for (j = 0; j < total; j++)     free(pts[j]);
    free(pts);
    
    return 1;
}


/*********************************************************************
* Author: Roummel F. Marcia
* Function: get_angles_slice
*
* This function computes the non-uniform angle distribution over half  
* an ellipse to determine where to place the ellipsoidal slices.
*
*********************************************************************/
void get_angles_slice(double thetas[], int nrings, double a, double b, double h)
{
    double  theta0, theta_targ, theta_start, arclength_start;
    double  hcirc, length_targ, length_start;
    int     k, n;

    theta0          = 0.0;
    theta_start     = theta0;
    arclength_start = 0.0;

    if (a >= b){
           if ( (nrings%2) == 0){
	          n = nrings/2;
		  length_start = h/2.0;
	   } else {
	          n = (nrings - 1)/2;
		  length_start = h;
	   }
    } else {
           length_start    = h;
           if ( (nrings%2) == 0)  n = (nrings/2) + 1;
	   else                   n = (nrings - 1)/2;
    }

    for (k = 0; k < n; k++)
    {
        if (a >= b){
	       length_targ = length_start + h*(double) k;
	       theta_targ  = compute_theta( a, b, theta_start, arclength_start, theta0, length_targ);
	       if (nrings%2 == 0)
	       {
		       thetas[n+k]   = (M_PI/2.0) + theta_targ;
		       thetas[n-k-1] = (M_PI/2.0) - theta_targ;
	       }
	       else
	       {
		       thetas[n+k+1] = (M_PI/2.0) + theta_targ;
		       thetas[n-k-1] = (M_PI/2.0) - theta_targ;
	       }
	       theta_start     = theta_targ;
	       arclength_start = length_targ;
	} else {
	       length_targ        = length_start + h*(double) k;
	       theta_targ         = compute_theta( b, a, theta_start, arclength_start, theta0, length_targ);
	       thetas[k]          = theta_targ;
	       thetas[nrings-k-1] = M_PI-theta_targ;
	       theta_start        = theta_targ;
	       arclength_start    = length_targ;
	}
    }

    if (nrings%2 == 1)  thetas[n] = M_PI/2.0;

    return;
}


/*********************************************************************
* Author: Roummel F. Marcia
* Function: get_angles_ellipse
*
* This function computes the non-uniform angle distribution over an 
* ellipse (with staggering).
*
*********************************************************************/
void get_angles_ellipse( double a, double b, double thetas[], int stagger, double nodes)
{
    int    i, n;
    double theta0, theta_start, length_targ, theta_targ;
    double arclength_start;
    double circ, hr, denom;
    hr    = pow(a-b, 2.0)/pow(a+b, 2.0);
    denom = 10.0 + sqrt( 4.0-(3.0*hr) );
    circ  = M_PI * (a+b) * ( 1.0+(3.0*hr/denom) );       // Ramanujan II for an ellipse
    theta0          = -M_PI/2.0;
    theta_start     =  theta0;
    arclength_start =  0.0;
    if (stagger == 1)
    {
        if ( (int) nodes%2 == 0 )  n =  ((int) nodes)/2;
        else                       n = (((int) nodes)-1)/2;
        for (i = 0; i < n; i++)
        {
            length_targ     = ( (double) i + 0.5 ) * circ / nodes;
            theta_targ      = compute_theta( a, b, theta_start, arclength_start, theta0, length_targ);
            thetas[i]       = theta_targ - theta0;
            theta_start     = theta_targ;
            arclength_start = length_targ;
        }
        if ( (int) nodes%2 == 0)
        {
            for (i = 0; i < n; i++) thetas[n+i]   = (2.0*M_PI) - thetas[n-(i+1)];
        }
        else
        {
            thetas[n] = M_PI;
            for (i = 0; i < n; i++) thetas[n+i+1] = (2.0*M_PI) - thetas[n-(i+1)];
        }
    }
    else
    {
        if ( (int) nodes%2 == 0 )    n = ( (int) nodes )/2;
        else                         n = (((int) nodes) + 1)/2;
        thetas[0] = 0.0;
        for (i = 1; i < n; i++)
        {
            length_targ     = ( (double) i ) * circ / nodes;
            theta_targ      = compute_theta ( a, b, theta_start, arclength_start, theta0, length_targ);
            thetas[i]       = theta_targ - theta0;
            theta_start     = theta_targ;
            arclength_start = length_targ;
        }
        if ( (int) nodes%2 == 0)
        {
            thetas[n] = M_PI;
            for (i = 1; i < n; i++)      thetas[n+i] = (2.0*M_PI) - thetas[n-i];
        }
        else
        {
            for (i = 0; i < (n-1); i++)  thetas[n+i] = (2.0*M_PI) - thetas[n - (i+1)];
        }
    }
    return;
}
/*********************************************************************
* Author: Roummel F. Marcia
* Function: compute_theta
*
* Given an initial angle and a specified arclength, this function
* computes the corresponding change in theta.
*
*********************************************************************/
double compute_theta( double a, double b, double theta_start, double arclength_start, double theta0, double length_targ)
{
    int     n;
    double  theta1, theta_a, theta_line, theta_quad, al, al1;
    double  lt0, lt1;
    double  circ, hr, denom;
    hr    = pow(a-b, 2.0)/pow(a+b, 2.0);
    denom = 10.0 + sqrt( 4.0-(3.0*hr) );
    circ  = M_PI * (a+b) * ( 1.0+(3.0*hr/denom) );       // Ramanujan II for an ellipse
    if (length_targ > (circ/2.0))
    {
        fprintf(stderr,"The target length is too long.\n");
        exit(1);
    }
    if (theta_start < theta0)   theta_start = theta0;
    theta1 = theta_start;
    lt1    = arclength_start;
    n      = 0;
    if (arclength_start >= length_targ) fprintf(stderr,"Ooops, something's wrong\n");
    // find interval [theta_a, theta1] that theta_targ is in
    while (lt1 < length_targ && n <= 90)
    {
        theta_a = theta1;
        theta1  = theta1 + M_PI/90;                    // go by increments of 2 degrees
        lt0     = lt1;
        lt1     = arclength( a, b, theta0, theta1);
        n++;
    }
    if (fabs(length_targ - lt1) < 0.00000001)
    {
        theta_line = theta1;
        theta_quad = theta1;
        al1        = arclength( a, b, theta0, theta_quad);
        al         = al1;
    }
    else
    {
        theta_line = theta_a + (theta1 - theta_a)*(length_targ - lt0)/(lt1 - lt0);     // linear interpolation between theta_a and theta1
        al1        = arclength( a, b, theta0, theta_line);
        // quadratic interpolation
        if (fabs( al1 - length_targ) < 0.00000001)
        {
            theta_quad = theta_line;
            al         = al1;
        }
        else
        {
            theta_quad = interp_quad( theta_a, theta_line, theta1, lt0, al1, lt1, length_targ);
            al         = arclength( a, b, theta0, theta_quad);
        }
    }
    /*
    if (fabs(al - al1) > 0.00001)
    {
        printf("Iterations = %4d\n", n);
        printf("theta_start= %12.8lf\n", theta_start);
        printf("theta_a    = %12.8lf\n", theta_a);
        printf("theta_line = %12.8lf\n", theta_line);
        printf("theta_quad = %12.8lf\n", theta_quad);
        printf("theta1     = %12.8lf\n\n", theta1);
        printf("lt0        = %12.8lf\n", lt0);
        printf("al1        = %12.8lf\n", al1);
        printf("al         = %12.8lf\n", al);
        printf("al_targ    = %12.8lf\n", length_targ);
        printf("lt1        = %12.8lf\n\n", lt1);
    }
    */
    return theta_quad;
}
/*********************************************************************
* Author: Roummel F. Marcia
* Function: interp_quad
*
* Given three data points, this function computes an approximate x
* value corresponding to a specified function value using
* quadratic approximation.
*
*********************************************************************/
double interp_quad( double t1, double t2, double t3, double r1, double r2, double r3, double al)
{
    double a, b, c;
    double theta;
    double det, rt1, rt2;
    det = (t2 - t1)*(t3 - t1)*(t3 - t2);
    if (det == 0.0)
    {
        fprintf(stderr,"Determinant is zero!\n");
        exit(1);
    }
    b = ( ((t3*t3) - (t1*t1))*(r2-r1) - ((t2*t2)-(t1*t1))*(r3-r1) )/det;
    a = ( -(t3 - t1)*(r2-r1)          + (t2 - t1)*(r3-r1))/det;
    c = r1 - (t1*b) - (t1*t1*a);
    rt1 = (-b + sqrt(b*b - 4.0*a*(c-al)))/(2.0*a);
    rt2 = (-b - sqrt(b*b - 4.0*a*(c-al)))/(2.0*a);
    if (rt2 >= t1 && rt2 <= t3)        theta = rt2;
    else                               theta = rt1;

    return theta;
}
/*********************************************************************
* Author: Roummel F. Marcia
* Function: arclength
*
* Given the axes a > b and two angles (theta0 < theta1), this function 
* computes the corresponding ellipsoidal arclength.
*
*********************************************************************/
double arclength( double a, double b, double theta0, double theta1)
{
    double  al = 0.0;
    double  e  = sqrt(1.0 - ((b/a)*(b/a)));
    double  intsin, nk;
    int     n;
    int     Nmax = 50;
    n = 1;
    intsin = integrate_sin_2n(n, theta0, theta1);
    al     = (theta1 - theta0);         // first term
    al     = al - (e*e*intsin/2.0);     // second term
    for (n = 2; n < Nmax; n++)
    {
        intsin = integrate_sin_2n(n, theta0, theta1);
        nk     = n_choose_k(2*n, n);
        al     = al + nk*(pow( e*e/4.0, (double) n))*intsin/(1.0 - 2.0*(double) n);
    }
    al = al*a;
    return al;
}
/*********************************************************************
* Author: Roummel F. Marcia
* Function: integrate_sin_2n
*
* Given two angles and an exponent n, this function computes the
* integral of (sin x)^2n.
*
*********************************************************************/
double integrate_sin_2n(int n, double theta0, double theta1)
{
    double int_val;
    int    k;
    double lval, uval;
    double intsum;
    double nk;
    nk      = n_choose_k(2*n, n);
    lval    = theta0 * nk / (pow(2.0, 2.0*(double) n));
    intsum  = 0.0;
    for (k = 0; k < n; k++)
    {
        nk     = n_choose_k(2*n, k);
        intsum = intsum + ( pow( -1.0, (double) k + (double) n) * nk * sin( 2.0*( (double) n - (double) k)*theta0) / ( (double) n - (double) k ) );
    }
    intsum  = intsum/( pow( 2.0, 2.0*(double) n) );
    lval    = lval + intsum;
    nk      = n_choose_k(2*n, n);
    uval    = theta1*nk/(pow(2.0, 2.0*(double) n));
    intsum  = 0.0;
    for (k = 0; k < n; k++)
    {
        nk     = n_choose_k(2*n, k);
        intsum = intsum + ( pow( -1.0, (double) k + (double) n) * nk * sin( 2.0*( (double) n - (double) k)*theta1) / ( (double) n - (double) k ) );
    }
    intsum  = intsum/( pow( 2.0, 2.0*(double) n) );
    uval    = uval + intsum;
    int_val = uval - lval;
    return int_val;
}
/*********************************************************************
* Author: Roummel F. Marcia
* Function: n_choose_k
*
* This function computes n choose k.  Returns a double from int n
* and k.
*
*********************************************************************/
double n_choose_k(int n, int k)
{
    double nval = 1.0;
    int    i;
    for (i = 0; i < k; i++)
    {
        nval = nval / ( (double) (k-i) );    // divide then multiply to prevent overflow
        nval = nval * ( (double) (n-i) );
    }
    return nval;
}


/*********************************************************************
 * Author: Roummel F. Marcia
 * Function: order_abc
 *
 * This function orders a1, b1, and c1 to a >= b >= c.
 *
 *********************************************************************/
void order_abc(double x, double y, double z, double *a, double *b, double *c, int index[])
{
     if (x >= y && x >= z){
           *a = x;
	   index[0] = 0;
	   if (y >= z){
	         *b = y;
		 *c = z;
		 index[1] = 1;
		 index[2] = 2;
	   } else {
	         *b = z;
		 *c = y;
		 index[1] = 2;
		 index[2] = 1;
	   }
     } else {
           if (y >= z){
	         *a = y;
		 index[0] = 1;
		 if (x >= z){
		      *b = x;
		      *c = z;
		      index[1] = 0;
		      index[2] = 2;
		 } else {
		      *b = z;
		      *c = x;
		      index[1] = 2;
		      index[2] = 0;
		 }
	   } else {
	         *a = z;
		 index[0] = 2;
		 if (x >= y){
		      *b = x;
		      *c = y;
		      index[1] = 0;
		      index[2] = 1;
		 } else {
		      *b  = y;
		      *c = x;
		      index[1] = 1;
		      index[2] = 0;
		 }
	   }
           
     }
     
     return;
}
