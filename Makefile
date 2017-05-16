EXECUTABLE=linearwaves
LIBRARY=liblinearwaves.so
SOURCES=coeffs.c disk.c main.c output.c ctridiag.c fft.c read_params.c torques.c viscosity.c linearwaves.c interpolation.c second_order.c
HEADER=linearwaves.h structs.h prototypes.h defines.h
LIBHEADER=$(HEADER) liblinear.h

LIBDIR=lib/
INCDIR=include/

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


nopara: $(CSOURCES) $(EXECUTABLE)
lib: CFLAGS += -fPIC
lib: $(CSOURCES) $(LIBRARY) resources

para: CFLAGS += -D_MPI
para: CC=$(MPICC)
para: $(CSOURCES) $(PEXECUTABLE)

install: lib installresources

all: $(CSOURCES) $(EXECUTABLE)

$(PEXECUTABLE): $(COBJECTS)
	$(MPICC)  $(COBJECTS) $(LDLIB) $(LDFLAGS) -o $@

$(EXECUTABLE): $(COBJECTS)
	$(CC)  $(COBJECTS) $(LDLIB) $(LDFLAGS) -o $@

$(LIBRARY): $(COBJECTS)
	$(CC) -shared $(COBJECTS) $(LDLIB) $(LDFLAGS) -o $@

installresources: installdirectories
	@cp $(LIBDIR)$(LIBRARY) /usr/local/lib/
	@cp $(SRC)prototypes.h /usr/local/include/liblinear/
	@cp $(SRC)defines.h /usr/local/include/liblinear/
	@cp $(SRC)structs.h /usr/local/include/liblinear/
	@cp $(SRC)liblinear.h /usr/local/include/liblinear/

installdirectories:
	@mkdir -p /usr/local/include/liblinear/

resources: directories
	@mv $(LIBRARY) $(LIBDIR)
	@cp $(SRC)prototypes.h $(INCDIR)
	@cp $(SRC)defines.h $(INCDIR)
	@cp $(SRC)structs.h $(INCDIR)
	@cp $(SRC)liblinear.h $(INCDIR)

directories:
	@mkdir -p $(LIBDIR)
	@mkdir -p $(INCDIR)



$(BIN)%.o: $(SRC)%.c $(CHEADER)
	$(CC) $(INCLIB) $(CFLAGS) $< -o $@


clean:
	rm $(COBJECTS)  $(EXECUTABLE) $(LIBDIR)$(LIBRARY) $(INCDIR)prototypes.h $(INCDIR)structs.h $(INCDIR)liblinear.h
