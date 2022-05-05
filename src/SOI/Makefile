###############################################################################
# Platform-dependent parameters
CC	= gcc
CPP	= g++
F77	= f77
MPIC	=
MPIF	= hf77
CFLAGS	= -g -O2 
FFLAGS	= -O2
LDFLAGS	=  
LIBS	= -lm
SYS	=
SPECIAL	= 
RANLIB	= /bin/true
SHELL	= /bin/sh
###############################################################################

TARGETS	=	test_energy test_coverage test_separation sample_icosa sample_sphere sample_ellipsoid sample_soi pts2pdb normalize set_sizes

UTIL =	rotutils.h	

COV 	=	test_coverage.o rotutils.o 
ENG 	=	test_energy.o rotutils.o  
SEP 	=	test_separation.o rotutils.o

ICOSA 	=	sample_icosa.o rotutils.o 
ELLIPSE 	=	sample_ellipsoid.o rotutils.o  
SPHERE 	=	sample_sphere.o rotutils.o
ROTS 	=	sample_soi.o rotutils.o
PTS 	=	pts2pdb.o rotutils.o
NORM	=	normalize.o rotutils.o
SIZES	=	set_sizes.o rotutils.o

#: Default rule
default: install clean   

# Rules to link object files into executables

test_coverage: $(COV) $(UTIL)
	$(CC) -o $@ $(COV) $(LDFLAGS) $(LIBS)

test_energy: $(ENG) $(UTIL)
	$(CC) -o $@ $(ENG) $(LDFLAGS) $(LIBS)

test_separation: $(SEP) $(UTIL)
	$(CC) -o $@ $(SEP) $(LDFLAGS) $(LIBS)

sample_sphere: $(SPHERE) $(UTIL)
	$(CC) -o $@ $(SPHERE) $(LDFLAGS) $(LIBS)

sample_icosa: $(ICOSA) $(UTIL)
	$(CC) -o $@ $(ICOSA) $(LDFLAGS) $(LIBS)

sample_ellipsoid: $(ELLIPSE) $(UTIL)
	$(CC) -o $@ $(ELLIPSE) $(LDFLAGS) $(LIBS)

sample_soi: $(ROTS) $(UTIL)
	$(CC) -o $@ $(ROTS) $(LDFLAGS) $(LIBS)

pts2pdb: $(PTS) $(UTIL)
	$(CC) -o $@ $(PTS) $(LDFLAGS) $(LIBS)

normalize: $(NORM) $(UTIL)
	$(CC) -o $@ $(NORM) $(LDFLAGS) $(LIBS)

set_sizes: $(SIZES) $(UTIL)
	$(CC) -o $@ $(SIZES) $(LDFLAGS) $(LIBS)

# Rule to build object file.o from file.c
.c.o:
	$(CC) -c $(CFLAGS) $(MYINC) $<

# Rule to build object file.o from file.C
.C.o:
	$(CPP) -c $(FFLAGS) $(MYINC) $<

# Rule to build object file.o from file.f
.f.o:
	$(F77) -c $(FFLAGS) $<

install: $(TARGETS)
	for executable in $(TARGETS); do \
		/bin/mv $${executable} ../bin/$${executable}; \
		chmod 755 ../bin/$${executable}; \
	done

clean:
	/bin/rm -f *.o *~
	/bin/rm -f $(TARGETS) core 

realclean: clean
	for executable in $(TARGETS); do \
                /bin/rm -f ../bin/$${executable}*; \
	done
