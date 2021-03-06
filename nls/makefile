FC = gfortran
CC = gcc

ZIP = zip

TARGET = solve
TARGETMOD = nls
TARGETLIB = native
VERSION = 0.1.1

BINDIR = ../bin
OBJDIR = ../obj
TARGETDIR = $(BINDIR)
DIRS = $(BINDIR) $(OBJDIR)

CFLAGS = -O3 -Wall -J $(OBJDIR)
LFLAGS = -std=f90 -lblas

SOURCES = nls.f90 \
		  solve.f95
OBJECTS = $(SOURCES:.f95=.o)

vpath %.o ../obj

all: checkdirs $(TARGET)

$(TARGET): $(addprefix $(OBJDIR)/, $(OBJECTS))
	$(FC) -L../obj $(LFLAGS) -o $(TARGETDIR)/$(TARGET) -J $(OBJDIR) $^

checkdirs: $(DIRS)

$(OBJDIR):
	mkdir $@

$(BINDIR):
	mkdir $@

$(OBJDIR)/%.o: %.f90
	$(FC) -c $(CFLAGS) -o $(OBJDIR)/$@ $<

install:

lib: glue

glue: checkdirs
ifeq ($(FCOMPILER),intelem)
	f2py2 --quiet -lmkl_def -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -c nls.f90 -m $(TARGETLIB) --fcompiler=intelem --f90flags='-mkl -heap-arrays'
else
	f2py2 --quiet -lblas -c nls.f90 -m $(TARGETLIB)
endif

clean:
	rm -rfv $(patsubst %, $(OBJDIR)/%, $(OBJECTS) $(TARGETMOD).mod)
	rm -rfv $(addprefix $(TARGETDIR)/, $(TARGETLIB).so $(TARGET) $(TARGETDIR)/nls-$(VERSION).zip)

.PHONY: all checkdirs clean install glue
.SUFFIXES: .f90 .o
