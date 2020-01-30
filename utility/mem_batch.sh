#! /bin/bash

#
# mem_batch.sh
#
# accepts commands over stdin, groups and runs them in batch using sge
#
# Copyright 2018 Peter Andrews @ CSHL
#

name=${1:-batch}
n_per_job=${2:-10}
n_jobs_max=${3:-2000}

n_jobs_at_start=$(jobids.sh -C)
n_jobs_to_run=$((n_jobs_max - n_jobs_at_start))
n_jobs=0
n_in_job=0
commands="echo"
mem_res=$( [ $(hostname -s) = mumdel ] && echo h_vmem || echo vf)

dir=${name}_sge
mkdir -p $dir
max_mem=10
while read mem command ; do
    if [ $n_in_job -ge $n_per_job ] || [ "$mem" = DONE ] ; then
        [ $n_in_job = 0 ] && break
        jname=b$n_jobs.$name
        out=$dir/$jname.txt
        err=$dir/$jname.err.txt
        echo $n_jobs $n_jobs_to_run $n_in_job $n_jobs_at_start $n_jobs_max
        echo "$commands ; echo All Done with batch job $jname" |
        qsub -cwd -l cpuver=5 -o $out -e $err -N $jname -l $mem_res=${max_mem}G > /dev/null ||
        break
        commands="echo"
        n_in_job=0
        max_mem=10
        n_jobs=$((n_jobs+1))
        [ "$mem" = DONE ] && break
    fi
    [ $n_jobs = $n_jobs_to_run ] && break
    commands="$command ; $commands"
    if [ $max_mem -lt $mem ] ; then
        max_mem=$mem
    fi
    n_in_job=$((n_in_job + 1))
done < <(cat ; echo DONE)

