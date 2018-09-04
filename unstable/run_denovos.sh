#! /bin/bash

. ~/mumdex/utility/shell.sh

ref=$1
families=$2
bridges=$3
samples=$4
chr=$5
start=$6
stop=$7
n_threads=$8
if [ "$n_threads" = "" ] || [ $n_threads -lt 1 ] ; then
    error usage run_denovos.sh ref families bridges samples chr start stop n_threads
fi

node=$(hostname -s)
echo running MUMdex population_denovos on $node

start_time=$(date +%s)

[ ! -e cand.txt ] || error output file cand.txt already exists

echo rsync
rsync -a ~/mumdex/ mumdex-code > /dev/null

echo compile
(cd mumdex-code ; make clean ; make -j 4 population_denovos) > /dev/null ||
error problem compiling MUMdex

compile=$(date +%s)

export PATH=mumdex-code:$PATH
set -o pipefail

echo
echo running population_denovos
echo population_denovos $ref $families $bridges $samples $n_threads $chr $start $stop
if population_denovos $ref $families $bridges $samples $n_threads $chr $start $stop > cand.txt ; then
    echo population_denovos succeeded
else
    errors population_denovos failed
fi

denovos=$(date +%s)

echo
(cd mumdex-code ; make clean > /dev/null)

clean=$(date +%s)

echo timing compile $((compile - start_time)) denovos $((denovos - compile)) clean $((clean-denovos)) all $((clean-start_time))

echo all done with population_denovos

exit 0

