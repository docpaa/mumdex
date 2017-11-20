#! /bin/bash

code_source=${1%/}/
code_target=${2%/}

shift 2

echo compiling programs $@ of $code_source in $code_target

(cd $code_source ; make clean)

rsync -aL $code_source $code_target 

(cd $code_target ; make -j $@ > /dev/null) && echo compiled



