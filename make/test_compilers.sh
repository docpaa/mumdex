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
out_dir=~/.compile-$(basename $PWD)
code_dir=$out_dir/code
rsync -aL --delete $PWD/ $code_dir
(cd $code_dir ; make clean > /dev/null)
mkdir -p $out_dir

function compile() {
    (
        compiler=$1
        nproc=$2
        export CXX=g++
        if [ $compiler = clang ] ; then
            export CXX=clang++
        elif [ $compiler = mac_clang ] ; then
            remote=andmin
            export CXX=clang++
        elif [ $compiler = 4.9.2 ] ; then
            remote=drpa.us
        elif [ $compiler = cygwin ] ; then
            remote=andrewsp@hannah5
        elif [ $compiler = default ] ; then
            a=b
        else
            export COMPILER_DIR=/data/software/gcc/$compiler
            export PATH=$COMPILER_DIR/bin:$PATH
        fi
        cd $code_dir
        if [ -z "$one_compiler" ] || [ $compiler = "$one_compiler" ] ; then
            start=$(date +%s)
            out="$(printf "Compiler %-9s" $compiler)"
            if [ ! -z $remote ] ; then
                rsync -aL --delete $PWD/ $remote:test_compilers
                ssh $remote "(hostname -s ; $CXX --version | head -n 1) 1>&2 && cd test_compilers && make clean > /dev/null && make -j 12 $program && ./even /dev/null" > /dev/null 2> $out_dir/$compiler.txt
            else
                echo -n "$out"
                ((hostname -s ; $CXX --version | head -n 1) 1>&2 && make -j $nproc $program && ./even /dev/null && make clean) > /dev/null 2> $out_dir/$compiler.txt
            fi
            sss=$?
            stop=$(date +%s)
            echo "$([ ! -z $remote ] && echo -n "$out")" $((stop-start)) s $sss e $(tail -n +3 $out_dir/$compiler.txt | wc -l) l $(head -n 2 $out_dir/$compiler.txt) $out_dir/$compiler.txt
        fi
    )
}

compile default $n_threads
for compiler in 9.3.0 8.4.0 7.4.0 6.5.0 5.5.0 ; do
    compile $compiler $n_threads
done
compile clang $n_threads
compile 4.9.2 12 &
compile mac_clang $n_threads &
wait
# compile cygwin 12
