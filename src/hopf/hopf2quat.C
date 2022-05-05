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

bool hopf2quat(vector < vector <double> > Points, bool writeq)
{
	double x1=0,x2=0,x3=0,x4=0;
    ofstream output;
    if (writeq)
    {
        output.open("data.qua");
        output << (int)Points.size() << "\n";
        output.precision(16);
    }
	for(int i=0;i<Points.size();i++)
	{
		x4=sin(Points[i][0]/2)*sin(Points[i][1]+Points[i][2]/2);
		x1=cos(Points[i][0]/2)*cos(Points[i][2]/2);
		x2=cos(Points[i][0]/2)*sin(Points[i][2]/2);
		x3=sin(Points[i][0]/2)*cos(Points[i][1]+Points[i][2]/2);
		if (writeq)
        {
            output << x1 << "\t" << x2 << "\t" << x3 << "\t" << x4 << endl;
            //output << Points[i][0]*180/M_PI << " " << Points[i][1]*180/M_PI << " " << Points[i][2]*180/M_PI << endl;
        }
	}
    if (writeq) output.close();
	return true;
}
