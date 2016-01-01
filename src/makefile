FC = gfortran
CC = gcc

TARGET = solve
CFLAGS = -O3 -Wall -J./
LFLAGS = -std=f95 -lblas
SOURCES = nls.f95 \
		  solve.f95
OBJECTS = $(subst '.f95', '.o', $(SOURCES))

all:
	$(FC) -c $(CFLAGS) $(SOURCES)
	$(FC) $(LFLAGS) -o $(TARGET) $(OBJECTS)

install:


glue:
	f2py2 -c $(SOURCE) -m $(TARGET)

clean:
	rm -rf *.mod *.o *.so $(TARGET)

.SUFFIXES: .f95 
