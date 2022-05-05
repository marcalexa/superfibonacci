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

int main(int argc, char *argv[])
{
	long int num_points=0;
	vector < double > Point;
	vector < vector<double> > Points;
	vector < vector<int> > Sequence_base;
	vector <int> temp;
	ifstream input;
	long int base_grid=0,cur_point=0;
	long int point_healpix=0;
	double point_S1=0;
	double theta=0,phi=0,psi=0;
	int limit=0;
	char filename[50];

//    cout << "Enter number of points in the sequence: ";
//    cin >> num_points;
    num_points = std::atoi(argv[1]);
//    cout << "Enter the filename in which sequence of the base grid is there: ";
//    cin >> filename;
	
	input.open("seq.txt");
	Sequence_base.resize(0);
	temp.resize(2);
	while(!input.eof())
	{
		input >> temp[0] >> temp[1];
		Sequence_base.push_back(temp);
	}
	input.close();
	Sequence_base.pop_back();
	
	//first seventy two points are the base grid points;
    std::cerr << "Creating ";
    std::cout << num_points << " ";
    std::cerr << "points by the Hopf method." << std::endl;

    auto tp0 = std::chrono::steady_clock::now();

	Points.resize(0);
	if(num_points<72)
		limit=num_points;
	else
		limit=72;

	for(int i=0;i<limit;i++)
	{
		Point.resize(0);
		pix2ang_nest(1,Sequence_base[i][0],&theta,&phi);
		switch(Sequence_base[i][1]) //mapping index on S1 to its angle value
		{
			case 0:
				point_S1=30;
				break;
			case 1:
				point_S1=90;
				break;
			case 2: 
				point_S1=150;
				break;
			case 3:
				point_S1=210;
				break;
			case 4: 
				point_S1=270;
				break;
			case 5: 
				point_S1=330;
				break;
		}
		psi=point_S1*M_PI/180;
		Point.push_back(theta);
		Point.push_back(phi);
		Point.push_back(psi);
		Points.push_back(Point);
	}

	
	for(int i=0;i<num_points-72;i++) //this will only be called if points are more than 72.
	{
		Point.resize(0);
		base_grid=i%72;
		cur_point=i/72;
		point_healpix=4*Sequence_base[base_grid][0];
		switch(Sequence_base[base_grid][1]) //mapping index on S1 to its angle value
		{
			case 0:
				point_S1=30;
				break;
			case 1:
				point_S1=90;
				break;
			case 2: 
				point_S1=150;
				break;
			case 3:
				point_S1=210;
				break;
			case 4: 
				point_S1=270;
				break;
			case 5: 
				point_S1=330;
				break;
		}
		Point=find_point(Sequence_base[base_grid][0],cur_point,1,point_healpix,point_S1); //current point value,level,current point in healpix, current point for S1
		Points.push_back(Point);
	}
	
	if(hopf2quat(Points,(argc > 2)))
    {
        auto tp1 = std::chrono::steady_clock::now();
        std::cerr << "Done in ";
        std::cout << (std::chrono::duration <double, std::milli> (tp1-tp0).count());
        std::cerr << " ms";
        std::cout << std::endl;
		return 0;
    }
	else
	{
		cout << "Problem in converting to quaternions" << endl;
		return 0;
	}

}
	
