#! /bin/bash

delim=' \t'
if [ "$1" != "" ] ; then
    delim="$1"
fi

perl -pe 's/['"$delim"']+/\n/g'
