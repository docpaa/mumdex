#! /bin/bash

sample=$1
mumdex=~/analysis/mums/wg-output/samples/$sample/mumdex
node=$(hostname -s)
noden=${node#wigclust}

output_dir=/data/unsafe/paa/bridges/$sample/
if [ -e $output_dir ] ; then
    echo output dir $output_dir already exists 1>&2
    exit 1
fi
mkdir -p $output_dir
cd $output_dir

echo running $sample in /mnt/$node$output_dir

n_threads=12
/usr/bin/time /data/unsafe/paa/mums/mumdex/bridges $mumdex $n_threads

available=$(df -B 1000000000 . | tail -n 1 | awk '{print $4}')
if [ $available -lt 500 ] ; then
    cd ../
    most_space=$(df -B 1000000000 /mnt/wigclust*/data | grep -v wigclust23 | sort -k4n | tail -n 1 | awk '{print $6}' | perl -pe 's+^/data+/mnt/'$node'/data+')
    mv $sample $most_space/unsafe/paa/bridges/$sample
    ln -s $most_space/unsafe/paa/bridges/$sample ~/analysis/mums/bridges
else
    ln -s /mnt/$node$output_dir ~/analysis/mums/bridges
fi
