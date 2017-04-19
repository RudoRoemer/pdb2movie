#################################################
# Makefile for HBdilute
# February 17, 2005. BMH
##################################################

CPP = g++
HBDILUTE = hbdilute
CPPFLAGS = -Wall -ansi

all: $(HBDILUTE)

$(HBDILUTE): hbdilute.o postscript.o hbfunctions.o 
	$(CPP) -o $@ hbdilute.o postscript.o hbfunctions.o $(CPPFLAGS)

clean:
	\rm -f *.o *.d
