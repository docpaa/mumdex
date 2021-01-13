#! /bin/bash

for dir in ~/analysis/smash-copy-number/jason{1..5} ~/analysis/smash_runs/jason{6..31} ~/analysis/smash_runs/jason_mice ; do
    name=${dir##*/}
    if [ ! -e site/data/jason/$name ] ; then
        echo processing $dir at $(date)
        ./setup_cn.sh jason $dir
    fi
done
    
