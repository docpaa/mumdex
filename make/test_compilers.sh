#! /bin/bash

#
# test multiple compiler versions
# to see if there are any compilation warnings or errors
#

if [ ! -e make/test_compilers.sh ] ; then
    echo Should be run from a source distribution directory 1>&2
    exit 1
fi

one_compiler=$1

program=
unset COMPILER_DIR
n_threads=$(nproc)
code_name=$(basename $PWD)
out_dir=~/.compile-test
mkdir -p $out_dir
code_dir=$out_dir/$code_name
rsync -aL --delete $PWD/ $code_dir
(cd $code_dir ; make clean > /dev/null)

function compile() {
    (
        compiler=$1
        nproc=$n_threads
        if [[ "$compiler" == *clang* ]] ; then
            export CXX=clang++
        else
            export CXX=g++
        fi
        if [ $compiler = clang ] || [ $compiler = default ] ; then
            a=b
        elif [ $compiler = mac_clang ] ; then
            remote=andmin
            nproc=12
        elif [ $compiler = gcc_4.9.2 ] ; then
            remote=drpa.us
            nproc=12
        elif [ $compiler = cygwin ] ; then
            remote=andrewsp@hannah5
            nproc=6
        else
            export COMPILER_DIR=/data/software/${compiler/_/\/}
            export PATH=$COMPILER_DIR/bin:$PATH
        fi
        cd $code_dir
        if [ -z "$one_compiler" ] || [ $compiler = "$one_compiler" ] ; then
            start=$(date +%s)
            out="$(printf "Compiler %-12s" $compiler)"
            if [ ! -z $remote ] ; then
                rdir=.test_compiler_$code_name
                rsync -aL --delete $PWD/ $remote:$rdir
                ssh $remote "(hostname -s ; $CXX --version | head -n 1) 1>&2 && cd $rdir && make clean > /dev/null && make -j 12 $program && ./even /dev/null" > /dev/null 2> $compiler.txt
            else
                echo -n "$out"
                ((hostname -s ; $CXX --version | head -n 1) 1>&2 && make -j $nproc $program && ./even /dev/null && make clean) > /dev/null 2> $compiler.txt
            fi
            sss=$?
            stop=$(date +%s)
            echo -n "$([ ! -z $remote ] && echo -n "$out") "
            printf "%3s s %1s e %3s l %-40s %-8s %s\n" $((stop-start)) $sss "$(tail -n +3 $compiler.txt | wc -l)" "$code_dir/$compiler.txt" "$(echo $(head -n 1 $compiler.txt))" "$(echo $(head -n 2 $compiler.txt | tail -n 1))"
        fi
    )
}

for compiler in \
    clang_{4..13}.x \
          gcc_{4.9.2,5.5.0,6.5.0,7.5.0,8.4.0,9.3.0,10.3.0,11.1.0} \
          default clang mac_clang ; do
    compile $compiler
done
# compile cygwin
