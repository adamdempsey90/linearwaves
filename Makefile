EXECUTABLE=linearwaves
SOURCES=coeffs.c disk.c main.c output.c ctridiag.c fft.c init.c torques.c viscosity.c linearwaves.c interpolation.c
HEADER=linearwaves.h

LAPACKLIB=-llapack -lblas
OMPLIB=-lgomp
MATHLIB=-lm
GSLLIB=-lgsl
FFTWLIB=-lfftw3

LDFLAGS= $(MATHLIB) $(LAPACKLIB) $(FFTWLIB) $(GSLLIB)

CFLAGS=-c -Wall -O3 

INCLIB=
LDLIB=


BIN=bin/
SRC=src/


UNAME := $(shell echo $(USER))


CC=gcc
#INCLIB=-I/usr/local/include/
#LDLIB=-L/usr/local/lib/
INCLIB=
LDLIB=

ifeq ($(UNAME),jupiter)
CC=gcc-6
endif

CC=mpicc
#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))
CHEADER=$(addprefix $(SRC),$(HEADER))




all: $(CSOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(COBJECTS)
	$(CC)  $(COBJECTS) $(LDLIB) $(LDFLAGS) -o $@

$(BIN)%.o: $(SRC)%.c $(CHEADER)
	$(CC) $(INCLIB) $(CFLAGS) $< -o $@


clean:
	rm $(COBJECTS) $(EXECUTABLE)
