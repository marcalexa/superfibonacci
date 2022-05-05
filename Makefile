###############################################################################
# Platform-dependent parameters
CPP	= g++
GPP = /usr/local/bin/g++-11
CFLAGS	= -O3 --std=c++17
EIGEN = -I/usr/local/include/eigen3
MYINC = -I/usr/local/include
LTBB = -ltbb
LDPATH = -L/usr/local/lib
GMP = -lgmp
SRC = ./src
BIN = ./bin
###############################################################################

TARGETS	= s3di s3soi s3sample s3cs s3ct s3pc

default: $(TARGETS)

s3di: $(SRC)/s3discrepancy.cc
	$(GPP) -o $(BIN)/$@ $< $(CFLAGS) $(EIGEN) $(MYINC) $(LDPATH) $(LTBB)

s3diall: $(SRC)/s3diall.cc
	$(CPP) -o $(BIN)/$@ $< $(CFLAGS) $(EIGEN)

s3pc: $(SRC)/s3pc.cc
	$(GPP) -o $(BIN)/$@ $< $(CFLAGS) $(EIGEN) $(MYINC) $(LDPATH) $(LTBB)

s3sample: $(SRC)/s3sample.cc
	$(CPP) -o $(BIN)/$@ $< $(CFLAGS) $(EIGEN)

s3cs: $(SRC)/s3cs.cc
	$(CPP) -o $(BIN)/$@ $< $(CFLAGS) $(EIGEN) $(GMP)

s3ct: $(SRC)/s3ct.cc
	$(CPP) -o $(BIN)/$@ $< $(CFLAGS) $(EIGEN) $(GMP)

s3soi: $(SRC)/s3soi.cc $(BIN)/rotutils.o
	$(CPP) -o $@ $^ $(CFLAGS) $(EIGEN)

$(BIN)/rotutils.o: $(SRC)/SOI/rotutils.c
	$(CPP) -c $< -O3
	
clean:
	/bin/rm -f *.o *~
	/bin/rm -f $(addprefix $(BIN)/, $(TARGETS)) core 

