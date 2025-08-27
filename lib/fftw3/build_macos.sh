#!/bin/bash

# macOS-specific FFTW3 build script for MURaM
# This script builds FFTW3 with MPI support for both single and double precision

# Setup Environment
SCRIPT_DIR=$(pwd)
OUTPUT_DIR=$SCRIPT_DIR
C_COMPILER=mpicc  # Use MPI-enabled C compiler from Homebrew

echo "Building FFTW3 for macOS with MPI support..."
echo "Output directory: $OUTPUT_DIR"
echo "Compiler: $C_COMPILER"

# Download FFTW3 if not already present
if [ ! -f fftw-3.3.10.tar.gz ]; then
    echo "Downloading FFTW 3.3.10..."
    wget -nc ftp://fftw.org/pub/fftw/fftw-3.3.10.tar.gz
else
    echo "FFTW 3.3.10 archive already exists, skipping download"
fi

# Extract if directory doesn't exist
if [ ! -d fftw-3.3.10 ]; then
    echo "Extracting FFTW 3.3.10..."
    tar -xzf fftw-3.3.10.tar.gz
fi

cd fftw-3.3.10

# Clean any previous builds
echo "Cleaning previous builds..."
make clean 2>/dev/null

# Build single precision version with MPI
echo "Building single precision version..."
./configure CC=$C_COMPILER CFLAGS="-O3" \
    --prefix=$OUTPUT_DIR \
    --enable-mpi \
    --enable-threads \
    --enable-float \
    --enable-shared=no

make -j$(sysctl -n hw.ncpu)
make install

# Clean before building double precision
make clean

# Build double precision version with MPI
echo "Building double precision version..."
./configure CC=$C_COMPILER CFLAGS="-O3" \
    --prefix=$OUTPUT_DIR \
    --enable-mpi \
    --enable-threads \
    --enable-shared=no

make -j$(sysctl -n hw.ncpu)
make install

cd ..

echo "FFTW3 build complete!"
echo "Libraries installed in: $OUTPUT_DIR/lib"
echo "Headers installed in: $OUTPUT_DIR/include"

# Verify installation
if [ -f "$OUTPUT_DIR/lib/libfftw3_mpi.a" ] && [ -f "$OUTPUT_DIR/lib/libfftw3f_mpi.a" ]; then
    echo "✅ Build successful! Both single and double precision MPI libraries are present."
else
    echo "⚠️  Warning: Some expected libraries may be missing. Check the build output for errors."
fi