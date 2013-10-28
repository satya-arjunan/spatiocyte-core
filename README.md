spatiocyte-core
===============

Multithreaded Standalone Spatiocyte

Compiling:
g++ -O3 -march=corei7  spatiocyte.cpp Model.cpp Stepper.cpp Compartment.cpp Species.cpp Diffuser.cpp VisualLogger.cpp -I. -I/home/satya/root/include -L/home/satya/root/lib -lRandom
