FC = gfortran
CC = gcc

TARGET = solve
TARGETLIB = nls

BINDIR = ../bin
OBJDIR = ../obj
TARGETDIR = $(BINDIR)
DIRS = $(BINDIR) $(OBJDIR)

CFLAGS = -O3 -Wall -J $(OBJDIR)
LFLAGS = -std=f95 -lblas

SOURCES = nls.f95 \
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

$(OBJDIR)/%.o: %.f95
	$(FC) -c $(CFLAGS) -o $(OBJDIR)/$@ $<

install:

glue: checkdirs
	f2py2 -lblas -c nls.f95 -m $(TARGETLIB)
	@mv -v {,$(TARGETDIR)/}$(TARGETLIB).so

clean:
	rm -rfv $(patsubst %, $(OBJDIR)/%, $(OBJECTS) $(TARGET).mod)
	rm -rfv $(addprefix $(TARGETDIR)/, $(TARGETLIB).so $(TARGET))

.PHONY: all checkdirs clean install glue
.SUFFIXES: .f95 .o
