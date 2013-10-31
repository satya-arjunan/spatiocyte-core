spatiocyte-core
===============

A multithreaded standalone spatiocyte package.


Compiling
~~~~~~~~~

1. Install Randomlib
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

2. Compile spatiocyte-core
    * $ make (or 'make -jx' with x the number of cores)


Running
~~~~~~~
You can spatiocyte-core by:
  * $ ./spatiocyte-core
