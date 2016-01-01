FC = gfortran
CC = gcc

TARGET = test_nls
CFLAGS = -I../src -O3 -Wall
LFLAGS = -std=f95 -lblas -L../src ../src/nls.o
SOURCES = test_nls.f95

all:
	$(FC) -o $(TARGET) $(CFLAGS) $(LFLAGS) $(SOURCES)

clean:
	rm -rf *.mod *.o *.so $(TARGET)

.SUFFIXES: .f95 
