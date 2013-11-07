VIS=\
		Visualizer

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
CXXFLAGS = -O3 -march=native
CXX = g++
GUILIBS = $(shell pkg-config --libs gtkmm-2.4 gtkglextmm-x11-1.2 libpng)
GUIFLAGS = $(shell pkg-config --cflags gtkmm-2.4 gtkglextmm-x11-1.2) -I.
CPPFLAGS = -DG_DISABLE_DEPRECATED -DGDK_PIXBUF_DISABLE_DEPRECATED -DPNG_SKIP_SETJMP_CHECK
OBJECTS=${SRC:=.o}
SPATIOCYTE_CORE = spatiocyte-core
VISUALIZER = visualizer


all: $(SPATIOCYTE_CORE) $(VISUALIZER)

$(SPATIOCYTE_CORE): $(OBJECTS)
		$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

$(VISUALIZER):
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(GUIFLAGS) -o $@ $(VIS).cpp $(GUILIBS)

%.o: %.cpp
		$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o $@ $<

clean:
		rm -f $(SPATIOCYTE_CORE) $(VISUALIZER) *.o
