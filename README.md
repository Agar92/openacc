# openacc
GPU OPENACC PROBLEM

Question Tuesday 27 August 2019 16.00:
The compile line:
cmake . -DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++ -DCMAKE_CXX_FLAGS="-acc -mcmodel=medium -ta=tesla:cc30 -Mnollvm -Mcuda=cuda10.1" -DCMAKE_CXX_STANDARD=17 -DACC=ON -DCUDA=ON

(In CMakeLists.txt #define OPENACC and #define CUDA are performed)

The compiler is PGI 19.4 pgc++.
