spatiocyte-core
===============

A multithreaded standalone spatiocyte package.

1. Get spatiocyte-core
    * $ git clone https://github.com/satya-arjunan/spatiocyte-core.git

2. Install Randomlib from http://sourceforge.net/projects/randomlib/files/distrib/
    * $ tar xzvf RandomLib-1.8.tar.gz 
    * $ cd RandomLib-1.8
    * $ mkdir BUILD
    * $ cd BUILD
    * $ CXXFLAGS="-O3 -march=corei7" CFLAGS="-O3 -march=corei7" cmake -DCMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX=$HOME/root ..
    * $ make
    * $ make test
    * $ make examples
    * $ cd examples
    * $ ./RandomTime
    * $ cd ..
    * $ make install

3. Compile spatiocyte-core
    * $ make (or 'make -jx' with x the number of available CPU cores)

4. Run spatiocyte-core
    * $ ./spatiocyte-core

5. Visualize logged data
    * $ ./visualizer
