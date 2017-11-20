#! /bin/bash

cd /home/paa/analysis/mums/nygc_out

cand=candidates_microsatellite.txt

grep40=$(awk '{printf "-e ^%s\n", $1}' ../wg-families.txt)
not40=$(awk '{print $1}' $cand | grep -v $grep40 | sort | uniq)

getnum='{if ($33 > 0) {++n ; tot += $33} ; atot += $33} END { print NR, n, tot/n}'
#getnum='{if ($33 >= 10) {++n ; tot += $33} ; atot += $33} END { print NR, n, tot/n}'
#getnum='{if ($33 > 0 && $33 < 10) {++n ; tot += $33} ; atot += $33} END { print NR, n, tot/n}'

echo $(grep $grep40 $cand | awk "$getnum") '<-' 40 families > tmp

n_tries=1000
n=0;
while [ $n -lt $n_tries ] ; do
    grepother40=$(randomize_order <(echo $not40 | perl -pe 's/ /\n/g') | head -n 40 | awk '{printf "-e ^%s\n", $1}')
    grep $grepother40 $cand | awk "$getnum" >> tmp
    n=$((n+1))
done

sort -k3n tmp > tmp2

(echo n_cand n_with_ms avg_ms
head -n 3 tmp2
echo
echo ... $((n_tries - 6)) skipped ...
echo
tail -n 4 tmp2
echo
cat -n tmp2 | grep famil | perl -pe 's/^ +//'|
awk '{print $2,$3,$4,$5,$6,$7,$1}') | even_columns ' '

