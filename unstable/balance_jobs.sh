#! /bin/bash

qstat -j '*SSC*map' > qstat2.txt

cat qstat2.txt | grep hard_queue | sort | uniq -c | sed 's/all.q@//' | sort -k1n > temp2.txt

cat temp2.txt

if [ -s temp2.txt ] ; then

    qstat -s p | grep ".m paa" | awk '{print $3}' > pend.txt
    
    # Move and link directories
    cat qstat2.txt | grep -e ^cwd -e ^hard_que | awk '{print $2}' | while read dir ; do read queue ; newnode=${queue#all.q@} ; newnode=${newnode%.cshl.edu};  oldnode=${dir#/mnt/} ; oldnode=${oldnode%%/*} ; alink=$(readlink $dir) ; sample=${dir/*\/} ; if [ "$alink" == "" ] && [ "$newnode" != "$oldnode" ] && grep $sample pend.txt > /dev/null ; then echo $sample $dir $queue $oldnode $newnode ; rm -f /mnt/$newnode/data/unsafe/paa/mums-new/samples/$sample ; mv $dir /mnt/$newnode/data/unsafe/paa/mums-new/samples/ ; ln -s /mnt/$newnode/data/unsafe/paa/mums-new/samples/$sample $dir ; rm /data/safe/paa/analysis/mums-new/output/samples/$sample ; ln -s /mnt/$newnode/data/unsafe/paa/mums-new/samples/$sample  /data/safe/paa/analysis/mums-new/output/samples/ ; fi  ; done
    
    echo done moving and linking directories

    n_empty=$(head -n 1 temp2.txt | awk '{print $1}')
    empty=$(head -n 1 temp2.txt | awk '{print $3}')
    n_full=$(tail -n 1 temp2.txt | awk '{print $1}')
    full=$(tail -n 1 temp2.txt | awk '{print $3}')
    n_empty=$((n_empty+1))
    
    if [ "$n_full" -gt "$n_empty" ] ; then
        one_full=$(qstat -q all.q@$full | awk '{if ($5 == "qw") print $3}' | sed 's/.m//' | tail -n 1)
        if [ "$one_full" != "" ] ; then
            echo transfer $one_full from $full $n_full to $empty $n_empty
            qalter -q all.q@$empty $one_full.map
            qalter -q all.q@$empty $one_full.calc
        else
            echo Could not find a good host transfer
        fi
    else
        echo Jobs are well balanced
    fi
    
else
    echo No jobs to balance
fi

rm -f qstat2.txt temp2.txt pend.txt
