#!/bin/bash
set -ex

# 1. Clone and build libMesh
git clone https://github.com/libMesh/libmesh.git
cd libmesh
git submodule update --init --recursive
mkdir build
cd build
../configure --prefix=$PREFIX
make -j${CPU_COUNT}
make install
cd ../..

# 2. Build and install apec_m with CMake
mkdir -p build_apec_m
cd build_apec_m
cmake ../src -DCMAKE_PREFIX_PATH=$PREFIX -DCMAKE_INSTALL_PREFIX=$PREFIX
make -j${CPU_COUNT}
make install