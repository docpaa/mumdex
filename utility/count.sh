#! /bin/bash

if [ "$1" = "-r" ] ; then
    r=r
fi

sort |
uniq -c |
sort -k1n$r

