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
	f2py2 --quiet -lblas -c nls.f90 -m $(TARGETLIB)
	@cp -v $(TARGETLIB).so $(BINDIR)/

clean:
	rm -rfv $(patsubst %, $(OBJDIR)/%, $(OBJECTS) $(TARGETMOD).mod)
	rm -rfv $(addprefix $(TARGETDIR)/, $(TARGETLIB).so $(TARGET) $(TARGETDIR)/nls-$(VERSION).zip)

.PHONY: all checkdirs clean install glue
.SUFFIXES: .f90 .o
