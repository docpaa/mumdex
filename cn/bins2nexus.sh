#! /bin/bash

bins_file=$1
base=$(basename $bins_file)
base=${base%_bins_results.txt}
# echo base $base

ratios=$base.ratios.txt
segs=$base.segs.txt

(
    cat $bins_file |
        table_columns.pl chr start stop ratio seg_ratio |
        tail -n +2 |
        awk '{print $1,$2,$3,$4*2,$5*2}'
    echo STOP
) |
    (
        echo -e "chr\tstart\tstop\tratio"
        echo -e "chr\tchrpos_start\tchrpos_stop\tn_bins\tseg_cn" 1>&2
        last_cn=-1
        last_start=0
        last_stop=0
        last_chr=0
        n_bins=0
        while read chr start stop ratio seg_ratio ; do
            if [ $last_cn = -1 ] ; then
                last_chr=$chr
                last_start=$start
                last_cn=$seg_ratio
            fi
            if [ $chr = STOP ] || [ $seg_ratio != $last_cn ] ||
               [ $chr != $last_chr ] ; then
                echo -e "$last_chr\t$last_start\t$last_stop\t$n_bins\t$last_cn" 1>&2
                last_start=$start
                last_cn=$seg_ratio
                n_bins=0
            fi
            if [ $chr = STOP ] ; then
                break
            fi
            echo -e "$chr\t$start\t$stop\t$ratio"
            last_chr=$chr
            last_stop=$stop
            n_bins=$((n_bins + 1))
        done
    ) > $ratios 2> $segs

