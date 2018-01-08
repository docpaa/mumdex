#! /bin/bash

#
# easily test multiple compiler versions
# to see if there are any compilation warnings or errors
#

for compiler in 4.9.2/rtf 5.5.0 6.4.0 7.2.0 ; do
    echo testing compiler version $compiler
    make clean > /dev/null
    export GCC_DIR=/data/software/gcc/$compiler
    make -j 24 > /dev/null
done
echo retoring default compilation
unset GCC_DIR
make clean > /dev/null
make -j 24 > /dev/null
