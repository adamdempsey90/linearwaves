EXECUTABLE=linearwaves
LIBRARY=linearwaves.so
SOURCES=coeffs.c disk.c main.c output.c ctridiag.c fft.c read_params.c torques.c viscosity.c linearwaves.c interpolation.c
HEADER=linearwaves.h

LAPACKLIB=-llapack -lblas
OMPLIB=-lgomp
MATHLIB=-lm
GSLLIB=-lgsl
FFTWLIB=-lfftw3

LDFLAGS= $(MATHLIB) $(LAPACKLIB) $(FFTWLIB) $(GSLLIB)

CFLAGS=-c -O3 -g -Wall -Wextra 

PEXECUTABLE=$(EXECUTABLE)

INCLIB=
LDLIB=


BIN=bin/
SRC=src/


UNAME := $(shell echo $(USER))


CC=gcc
MPICC=mpicc


#INCLIB=-I/usr/local/include/
#LDLIB=-L/usr/local/lib/
INCLIB=
LDLIB=

ifeq ($(UNAME),jupiter)
CC=gcc-6
endif

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))
CHEADER=$(addprefix $(SRC),$(HEADER))


para: CFLAGS += -D_MPI
para: CC=$(MPICC)
para: $(CSOURCES) $(PEXECUTABLE)

nopara: $(CSOURCES) $(EXECUTABLE)

all: $(CSOURCES) $(EXECUTABLE)

$(PEXECUTABLE): $(COBJECTS)
	$(MPICC)  $(COBJECTS) $(LDLIB) $(LDFLAGS) -o $@

$(EXECUTABLE): $(COBJECTS)
	$(CC)  $(COBJECTS) $(LDLIB) $(LDFLAGS) -o $@

$(BIN)%.o: $(SRC)%.c $(CHEADER)
	$(CC) $(INCLIB) $(CFLAGS) $< -o $@


clean:
	rm $(COBJECTS) $(EXECUTABLE)
