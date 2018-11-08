#! /bin/bash

n_threads="$1"
if [ "$n_threads" = "" ] ; then
    n_threads=5
fi

bed=~/hg38/hg38.1M.bed
if [ "$2" != "" ] ; then
    bed="$2"
fi

out=out
if [ "$3" != "" ] ; then
    out="$3"
fi

n_jobs_max=2000 &&
n_jobs=$(jobids.sh C) &&
tac $bed |
while read chr start stop ; do
    dir=$out/$chr/$start-$stop
    name=$chr.$start
    if [ ! -e $dir ] ; then
        n_jobs=$((n_jobs + 1))
        if [ $n_jobs -le $n_jobs_max ] ; then
            mkdir -p $dir
            (
                cd $dir
                echo submitting job $n_jobs of $n_jobs_max in $dir &&
                qsub -N $name -l h_vmem=40G -pe smp $n_threads \
                    -cwd -o out.txt -e err.txt \
                    ~/mumdex/unstable/run_transmission.sh ~/hg38/hg38.fa \
                    ~/ssc_hg38/{families_2380.txt,transmission/bridges} \
                    $chr $start $stop 10 > /dev/null ||
                error qsub failed
            ) || exit 1
        else
            echo max queue length has been reached
            break
        fi
    fi
done

exit 0

tail -qn 1 out/*/*/out.txt | grep 'all done with transmission' | wc -l

tail -qf -n +2 out/chr{{1..22},X,Y}/*/cand.txt | even_columns ' ' yes | cat -n

while true ; do qstat | grep ' r '|tail -n 1 | ~/mumdex/utility/fixed.sh ; sleep 60 ; done | uniq


