#! /bin/bash

#
# easily test multiple compiler versions
# to see if there are any compilation warnings or errors
#
nproc=$(nproc)
mkdir -p ~/.compile

cd ~/mumdex
for compiler in 9.1.0 8.3.0 5.5.0 6.5.0 7.4.0 ; do
    echo testing compiler version $compiler
    make clean > /dev/null
    export GCC_DIR=/data/safe/paa/gcc/$compiler
    make -j $nproc > /dev/null 2> ~/.compile/$compiler.txt
    sss=$?
    cat ~/.compile/$compiler.txt
    echo $sss status and $(cat ~/.compile/$compiler.txt | wc -l) error/warning lines in ~/.compile/$compiler.txt
done

# special clang compile test
compiler=clang
cd ~/mumdex
echo testing compiler version $compiler
make clean > /dev/null
cd ../
rsync -av --progress mumdex/ andmin:mumdex > /dev/null
ssh andmin "cd mumdex && make clean > /dev/null && make -j 12" > /dev/null 2> ~/.compile/$compiler.txt
sss=$?
cat ~/.compile/$compiler.txt
echo $sss status and $(cat ~/.compile/$compiler.txt | wc -l) error/warning lines in ~/.compile/$compiler.txt

# special 4.9.2 compile test
compiler=4.9.2
cd ~/mumdex
echo testing compiler version $compiler
make clean > /dev/null
cd ../
rsync -av --progress mumdex/ drpa.us:mumdex > /dev/null
ssh drpa.us "export GCC_DIR=/data/software/gcc/$compiler/rtf && cd mumdex && make clean > /dev/null && make -j 12" > /dev/null 2> ~/.compile/$compiler.txt
sss=$?
cat ~/.compile/$compiler.txt
echo $sss status and $(cat ~/.compile/$compiler.txt | wc -l) error/warning lines in ~/.compile/$compiler.txt

exit 0
echo restoring default compilation
unset GCC_DIR
make clean > /dev/null
make -j $nproc > /dev/null
