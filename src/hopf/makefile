# -----------------------------------------------------------------------------
#
#  Copyright (C) 2009  Anna Yershova, Swati Jain, 
#                      Steven M. LaValle, Julie C. Mitchell
#
#
#  This file is part of the Incremental Successive Orthogonal Images (ISOI)
#
#  ISOI is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  ISOI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
#  For more information about ISOI see http://rotations.mitchell-lab.org/
#
#----------------------------------------------------------------------------- */

SO3_sequence: main.o find_point.o pix2ang_nest.o mk_pix2xy.o hopf2quat.o
	g++ -o SO3_sequence main.o find_point.o pix2ang_nest.o mk_pix2xy.o hopf2quat.o -O3
main.o: main.C so3_sequence.h
	g++ -c main.C -O3
find_point.o: find_point.C so3_sequence.h
	g++ -c find_point.C -lm -O3
hopf2quat.o: hopf2quat.C so3_sequence.h
	g++ -c hopf2quat.C -lm -O3
pix2ang_nest.o: pix2ang_nest.c
	gcc -c pix2ang_nest.c -lm -O3
mk_pix2xy.o: mk_pix2xy.c
	gcc -c mk_pix2xy.c -lm -O3

clean:
	rm -f *.o
	rm -f SO3_sequence
