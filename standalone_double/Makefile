CC=g++

CFLAGS= -O3 -g -std=c++0x #-Wall
LFLAGS= -L/usr/lib -lboost_program_options-mt

SRCFILES=cpu_driver.cpp constants.cpp io.cpp laser.cpp particle.cpp particle_list.cpp physobj.cpp double3.cpp arg_parser.cpp

OBJFILES=$(SRCFILES:.cpp=.o)
EXEFILE=double_driver

all: $(SRCFILES) $(EXEFILE)

$(EXEFILE):  $(OBJFILES)
	$(CC) $(LFLAGS) -o $@ $(OBJFILES)

.cpp.o:
	$(CC) -c $(CFLAGS) -o $@ $<

clean:
	rm -f *.o *~ $(EXEFILE)