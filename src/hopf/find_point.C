/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 2009  Anna Yershova, Swati Jain, 
 *                      Steven M. LaValle, Julie C. Mitchell
 *
 *
 *  This file is part of the Incremental Successive Orthogonal Images (ISOI)
 *
 *  ISOI is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  ISOI is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this software; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about ISOI see http://rotations.mitchell-lab.org/
 *
 *----------------------------------------------------------------------------- */

#include"so3_sequence.h"

vector < double >find_point(int base_grid, long int point,long int level,long int healpix_point,double s1_point)
{
	int position=point%8;
	long int quo=0;
	vector <double> Point;
	double interval=30/level;
	// the choosing of the order of the first resolution 4 points depends on which base healpix grid we are now dividing 

	if(base_grid == 6 or base_grid == 7)
	{
		switch(position) //this position tells which of the eight points of the cube to consider 
		{ 
			case 0:
				healpix_point+=3;
				s1_point-=(interval/2);
				break;
			case 1:
				healpix_point+=0;
				s1_point+=(interval/2);
				break;
			case 2: 
				healpix_point+=3;
				s1_point+=(interval/2);
				break;
			case 3:
				healpix_point+=0;
				s1_point-=(interval/2);
				break;
			case 4:
				healpix_point+=2;
				s1_point-=(interval/2);
				break;
			case 5:
				healpix_point+=1;
				s1_point+=(interval/2);
				break;
			case 6:
				healpix_point+=2;
				s1_point+=(interval/2);
				break;
			case 7:
				healpix_point+=1;
				s1_point-=(interval/2);
				break;
			
		}
	}
	else if(base_grid == 3 or base_grid == 1 or base_grid == 9 or base_grid == 11)
	{
		switch(position) 
		{
			case 0:
				healpix_point+=3;
				s1_point-=(interval/2);
				break;
			case 1:
				healpix_point+=0;
				s1_point+=(interval/2);
				break;
			case 2: 
				healpix_point+=3;
				s1_point+=(interval/2);
				break;
			case 3:
				healpix_point+=0;
				s1_point-=(interval/2);
				break;
			case 4:
				healpix_point+=1;
				s1_point-=(interval/2);
				break;
			case 5:
				healpix_point+=2;
				s1_point+=(interval/2);
				break;
			case 6:
				healpix_point+=1;
				s1_point+=(interval/2);
				break;
			case 7:
				healpix_point+=2;
				s1_point-=(interval/2);
				break;
			
		}
	}
	else if(base_grid == 2 or base_grid == 0 or base_grid == 10 or base_grid == 8)
	{
		switch(position) 
		{ 
			case 0:
				healpix_point+=0;
				s1_point-=(interval/2);
				break;
			case 1:
				healpix_point+=3;
				s1_point+=(interval/2);
				break;
			case 2: 
				healpix_point+=0;
				s1_point+=(interval/2);
				break;
			case 3:
				healpix_point+=3;
				s1_point-=(interval/2);
				break;
			case 4:
				healpix_point+=1;
				s1_point-=(interval/2);
				break;
			case 5:
				healpix_point+=2;
				s1_point+=(interval/2);
				break;
			case 6:
				healpix_point+=1;
				s1_point+=(interval/2);
				break;
			case 7:
				healpix_point+=2;
				s1_point-=(interval/2);
				break;
			
		}
	}
	else if(base_grid == 4 or base_grid == 5)
	{
		switch(position) 
		{ 
			case 0:
				healpix_point+=0;
				s1_point-=(interval/2);
				break;
			case 1:
				healpix_point+=3;
				s1_point+=(interval/2);
				break;
			case 2: 
				healpix_point+=0;
				s1_point+=(interval/2);
				break;
			case 3:
				healpix_point+=3;
				s1_point-=(interval/2);
				break;
			case 4:
				healpix_point+=2;
				s1_point-=(interval/2);
				break;
			case 5:
				healpix_point+=1;
				s1_point+=(interval/2);
				break;
			case 6:
				healpix_point+=2;
				s1_point+=(interval/2);
				break;
			case 7:
				healpix_point+=1;
				s1_point-=(interval/2);
				break;
			
		}
	}
		
	quo=point/8;
	if(quo==0)
	{
		long int nside=pow(2,level);
		double theta=0,phi=0,psi=0;
		pix2ang_nest(nside,healpix_point,&theta,&phi);
		psi=s1_point*M_PI/180;
		Point.resize(0);
		Point.push_back(theta);	
		Point.push_back(phi);	
		Point.push_back(psi);
		return Point;
	}
	else
	{
		return find_point(base_grid,quo-1,level+1,4*healpix_point,s1_point); 
	}
}	
		
