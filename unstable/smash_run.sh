#! /bin/bash


if [ $# != 6 ] ; then
    echo usage: smash_run.sh mumdex_dir ref_file bins_dir n_bins,... data_dir out_dir 1>&2
    exit 1
fi

mumdex_dir="$1"
ref="$2"
bins="$3"
n_bins="$4"
data_dir="$5"
out_dir="$5"

cd $out_dir

find $data_dir/ -name '*.R1.fastq.gz'

