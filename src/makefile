FC = gfortran
CC = gcc

TARGET = solve
TARGETLIB = nls

OBJECTSDIR = ../obj
CFLAGS = -O3 -Wall -J $(OBJECTSDIR)
LFLAGS = -std=f95 -lblas

SOURCES = nls.f95 \
		  solve.f95
OBJECTS = $(SOURCES:.f95=.o)

vpath %.o ../obj

all: checkdir $(TARGET)

$(TARGET): $(addprefix $(OBJECTSDIR)/, $(OBJECTS))
	$(FC) -L../obj $(LFLAGS) -o $(TARGET) -J $(OBJECTSDIR) $^

checkdir: $(OBJECTSDIR)

$(OBJECTSDIR):
	mkdir $@

$(OBJECTSDIR)/%.o: %.f95
	$(FC) -c $(CFLAGS) -o $(OBJECTSDIR)/$@ $<

install:


glue:
	f2py2 -lblas -c nls.f95 -m $(TARGETLIB)

clean:
	rm -rfv $(patsubst %, $(OBJECTSDIR)/%, $(OBJECTS))
	rm -rfv $(TARGETLIB).so $(TARGET)

.PHONY: all clean
.SUFFIXES: .f95 .o
