#! /bin/bash

#
# run bridge_figure script for all candidates
#

process_memory=30G

# machine-dependent setup
if [ "$SGE_ROOT" = /opt/sge6-2 ] ; then
    # wigclust
    options="-q all.q -l vf=$process_memory"
    base=~/analysis/mums
else
    # NYGC
    options="-l h_vmem=$process_memory"
    base=~
fi

code=~/mumdex
output=$base/bridges_out/figures
families=$base/wg-families.txt
samples=$base/samples
bridges=$base/bridges

# subjobs per job
max_subjobs_per_job=1

# make output directory
mkdir -p $output
cd $output

# Jobs directory
jobs=$output/jobs/job.$$
mkdir -p $jobs
subjobs=$jobs/to_run.txt

# Code for all runs
jobscode=$jobs/code
rsync -aL $code/ $jobscode
(cd $jobscode ; make clean > /dev/null)

# Select subjobs to run
cat $output/../*/*-*/out.txt | grep -v family |
sort -k1,1 -k7,7n -k9,9n | uniq |
awk '{print $8,$9,$10,$11,$12,$13,$16,$1,$3,$17,$14}' |
while read chr1 pos1 high1 chr2 pos2 high2 inv fam mem off type ; do
    thigh1=${high1/1/high}
    thigh1=${thigh1/0/low}
    thigh2=${high2/1/high}
    thigh2=${thigh2/0/low}
    out=$fam.$chr1.$pos1.$thigh1.$chr2.$pos2.$thigh2.${inv}inv.${off}off.$type
    pdf=$output/$out.pdf
    if [ ! -e $pdf ] ; then
        echo "bridge_figure $samples $families $chr1 $pos1 $high1 $chr2 $pos2 $high2 $inv $fam $mem $off $type 2> $out.err.txt > $out.out.txt"
    fi
done > $subjobs

if [ ! -e $subjobs ] ; then
    echo no subjobs found to run
    exit 1
fi

# split up subjobs
split -l $max_subjobs_per_job -d -a 5 $subjobs $jobs/job || exit 1
rm $subjobs

# set up sge and job scripts
for job in $jobs/job* ; do
    mv $job temp
    mkdir $job
    mv temp $job

    jobcode=$job/code

    jobscript=$job/job.sh
    cat $job/temp |
    while read subjob ; do
        echo $subjob >> $jobscript
    done
    chmod a+x $jobscript
    rm $job/temp

    sgescript=$job/sge.sh
    name=${job##*/}
    touch $job/sge.err.txt
    touch $job/sge.out.txt
    touch $job/status.txt
    cat <<EOF > $sgescript
open=\$(ulimit -a| grep open | awk '{print \$4}')
start=\$(date +%s)
echo starting bridge_figure job on \$(hostname) at \$(date) >> $job/status.txt
rsync -a $jobscode/ $jobcode
(cd $jobcode ; make bridge_figure > /dev/null) || exit 2
export PATH=$jobcode:$PATH
cd $output
$jobscript
status=\$?
echo finished bridge_figure job. status: \$status >> $job/status.txt
stop=\$(date +%s)
echo total time was \$((stop - start)) >> $job/status.txt
exit \$status
EOF
    chmod a+x $sgescript

    echo qsub $options -N $name -o $job/sge.out.txt -e $job/sge.err.txt $sgescript | tee -a $jobs/sge.sh
done
if false ; then
    bash $jobs/sge.sh
else 
    echo To submit jobs, run:
    echo bash $jobs/sge.sh
fi
echo jobdir is $jobs
