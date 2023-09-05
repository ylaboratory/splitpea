#!/bin/bash

wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
tar xf htslib-1.17.tar.bz2
rm -f htslib-1.17.tar.bz2
cd htslib-1.17
./configure --prefix=$HOME
make
make install
