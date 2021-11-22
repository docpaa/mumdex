#! /bin/bash

echo Should be 55 files per sample produced:
for dir in $(ls -d */ | grep -v data); do echo $dir $(ls $dir | wc -l) ; done | sort -k2nr

echo Should show no bad segments files:
for bin in $(echo $bins | sed 's/,/ /g'); do for dir in $(ls -d */ | grep -v data) ; do dir=${dir%/} ; if [ ! -e $dir/${dir}_${bin}_bins_segments.txt ] ; then echo $dir/${dir}_${bin}_bins_segments.txt ; fi ; if [ ! -e $dir/${dir}_${bin}_bins_segments.err.txt ] ; then echo $dir/${dir}_${bin}_bins_segments.err.txt ; fi ; done ; done

echo Should show no bad lines in segments error file:
for file in */*segments.err.txt ; do if ! grep done $file > /dev/null ; then echo $file ; fi ; done

echo Generating overall stats file:
(cat */*stats.txt | head -n 1 ; cat */*stats.txt | grep -v sample | sort | uniq ) | tee stats.txt 
