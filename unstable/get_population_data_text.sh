#! /bin/bash

bridges_dir=$1
ref=$2
chr=$3
start=$4
stop=$5

index=$(perl -ne 'print $n++ if /^'$chr'\s/' $ref.fai)

echo Info is $ref $bridges_dir $chr $start $stop

for sample in $(cd $bridges_dir/ ; ls -d SSC*/ | perl -pe 's+/++') ; do
    bridge_dir=$bridges_dir/$sample
    bridge_file=$bridge_dir/chrbridges.$index.bin
    window=$chr:$start-$stop
    mkdir -p $window
    out=$window/$sample.txt.gz
    ~/mumdex/bridges2txt $ref $bridge_file $chr $start $stop | gzip -c > $out
    echo $sample $(ls -lh $out)
done





