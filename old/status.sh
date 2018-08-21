#! /bin/bash
exit 0

for n in {1..5} ; do echo ; done

done=0
running=0
bad=0
waiting=0
for status in $(ls -drt */*/status.txt) ; do
    mod_at="$(s=$(grep time $status| awk '{print $4}') ; if [ "$s" = "" ] ; then echo x ; else echo $s ; fi) $(stat -c %y $status)"
    dir=${status%/status.txt}
    lines=$(cat $status | wc -l)
    if [ $lines -gt 0 ] ; then
	status=$(grep finished $status | awk '{print $4}')
	if [ "$status" = 0 ] ; then
	    err=$dir/err.txt
	    if grep finished $err > /dev/null ; then
		if grep 'all done' $err > /dev/null ; then 
       		    echo $dir ok $mod_at
		    done=$((done+1))
		else
		    echo $dir good status and finished but not done $mod_at
		fi
	    else
		echo $dir good status but not finished $mod_at
	    fi
	elif [ "$status" = "" ] ; then	    
	    echo $dir processing $mod_at
	    # ls -lrt $dir/out.txt
	    running=$((running+1))
	else
	    echo $dir bad $status $mod_at
	    bad=$((bad+1))
	fi
    else
	waiting=$((waiting+1))
    fi
done > tmp

cat tmp | ~/mumdex/even_columns ' ' | cat -n
echo 
qstat | grep $USER > tmp
finished=$(grep finished jobs/*/*/status.txt | wc -l)
echo jobs: $(grep -v qw tmp | wc -l) running $(grep qw tmp | wc -l) waiting $finished finished
echo processes: $done done $running running $bad bad $waiting waiting
if [ $finished = 0 ] ; then
    exit 0
fi
echo timing: $(grep time */*/status.txt | grep -v GL | awk '{print $4}' | sort -k1nr | awk '{tot+=$1; print} END {print "avg of ", NR, "is", tot/NR}' )

find . -name out.txt | xargs cat | grep -v family | awk '{print $3}' | sort | uniq -c | awk '{print ; tot += $1} END{print tot, "total"}' | sort -k1nr > tmp
pro=$(cat tmp | grep proband | awk '{print $1}')
sib=$(cat tmp | grep sibling | awk '{print $1}')
echo candidates: $(cat tmp) $(~/mumdex/bin_p $((pro + sib)) $pro)

find . -name out.txt | xargs cat | grep -v -e family | awk '{if ($4 <= 1 && $18 >= 5 && $27 >= 10 && $19 >= 25 && $20 >= 25 && $22 >= 20 && $24 >= 20) print }' > filt
cat filt | awk '{print $3}' | sort | uniq -c | awk '{print ; tot += $1} END{print tot, "total"}' | sort -k1nr > tmp
pro=$(cat tmp | grep proband | awk '{print $1}')
sib=$(cat tmp | grep sibling | awk '{print $1}')
echo filtered: $(cat tmp) $(~/mumdex/bin_p $((pro + sib)) $pro)

cat filt | awk '{print $14}' | grep -v type | sort | uniq -c | sort -k1n

echo SGE messages:
find jobs -name sge*.txt | xargs cat | sort | uniq -c | fgrep -v -e '***' -e unlimited -e localhost -e LEGAL -e bytes -e accessed -e prohibited -e UNAUTH -e 'core file' -e enforcement -e property -e approvals -e information -e permitted -e communication -e processes -e 'open files' -e 'L E G A L' -e 'pending signals' -e 'scheduling' -e 'time priority'
echo

(find . -name out.txt | xargs cat | head -n 1 ; cat filt) > filtered_candidates.txt

rm tmp filt
