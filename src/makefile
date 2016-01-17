FC = gfortran
CC = gcc

TARGET = solve
TARGETLIB = nls

OBJECTDIR = ../obj
CFLAGS = -O3 -Wall -J $(OBJECTDIR)
LFLAGS = -std=f95 -lblas


SOURCES = nls.f95 \
		  solve.f95
OBJECTS = $(SOURCES:.f95=.o)

all: checkdir
	$(FC) -c $(CFLAGS) $(SOURCES)
	$(FC) $(LFLAGS) -o $(TARGET) -J $(OBJECTDIR) $(OBJECTS)

checkdir: $(OBJECTDIR)

$(OBJECTDIR):
	mkdir $@

install:


glue:
	f2py2 -lblas -c nls.f95 -m $(TARGETLIB)

clean:
	rm -rf *.mod $(OBJECTS) *.so $(TARGET)

.SUFFIXES: .f95 
