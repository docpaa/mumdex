#! /bin/bash

bridges_dir=$1
ref=$2
chr=$3
start=$4
stop=$5

index=$(perl -ne '$n++ ; print $n if /^'$chr'\s/' $ref.fai)

echo Info is $ref $bridges_dir $chr $start $stop

window=$chr-$start-$stop
mkdir -p $window
cd $window

rsync -a --exclude .git ~/mumdex/ mumdex
(cd mumdex && make -j 4 bridges2txt)

n=0
for sample in $(cd $bridges_dir/ ; ls -d $(echo SSC*/) | perl -pe 's+/++g') ; do
    bridge_dir=$bridges_dir/$sample
    bridge_file=$bridge_dir/chrbridges.$index.bin
    out=$sample.txt.gz
    n=$((n+1))
    echo $window $n 1>&2
    if [ ! -e $out ] ; then
        echo "./mumdex/bridges2txt $ref $bridge_file $chr $start $stop | gzip -c > $out && echo finished $sample $n"
    fi
done | qsub -cwd -N test -o out.txt -e err.txt 
