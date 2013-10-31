SRC=\
		spatiocyte\
		Model\
		Stepper\
		Compartment\
		Species\
		Diffuser\
		VisualLogger

IFLAGS = -I. -I$(HOME)/root/include
LDFLAGS = -L$(HOME)/root/lib -lRandom
CXXFLAGS = -O3 -march=corei7
CXX = g++
OBJECTS=${SRC:=.o}
SPATIOCYTE_CORE = spatiocyte-core

all: $(SPATIOCYTE_CORE)

$(SPATIOCYTE_CORE): $(OBJECTS)
		$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)
		rm *.o

%.o: %.cpp
		$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o $@ $<

clean:
		rm -f $(SPATIOCYTE_CORE)
