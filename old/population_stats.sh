#! /bin/bash

min_bridge_count=5
min_ref=10
max_families=1
min_support=25
min_mate_support=20
max_population_people=1000
max_parent_coverage=2000000000

filter="s/,/ /g; s/bridges/mb fb pb sb/; s/anchorsA/maA faA paA saA/; s/coveragesA/mcA fcA pcA scA/; s/anchorsB/maB faB paB saB/; s/coveragesB/mcB fcB pcB scB/;"

# cat chr*/*/out.txt | sort | uniq | perl -pe "$filter" | tac | head -n 2 | ~/code/scripts/number_columns.sh

cat chr*/*/out.txt | sort | uniq | perl -pe "$filter" | grep -v family |
sort -k7,7n -k9,9n > all_cand.txt
#echo candidates loose $(cat all_cand.txt | wc -l)

cat all_cand.txt | awk '{if ($18 >= '$min_bridge_count' && $27 >= '$min_ref' && $19 >= '$min_support' && $20 >= '$min_support' && $22 >= '$min_mate_support' && $24 >= '$min_mate_support' && $40 <='$max_parent_coverage' && $41 <= '$max_parent_coverage') print}'  > denovo_cand.txt
#echo candidates denovo $(cat denovo_cand.txt | wc -l)

echo How many families was denovo seen in:
cat denovo_cand.txt | awk '{print $8,$9,$10,$11,$12,$13,$16}' | sort | uniq -c | awk '{print $1}' | sort | uniq -c

header=$(cat chr*/*/out.txt | head -n 1 | perl -pe "$filter")

if true ; then
    (echo A B C CovChild BridgesChild BridgesNonFamily NNonFamily $header ;
        cat denovo_cand.txt | awk '{print ($4 == 1), ($54 == 0), ($48 == 0), ($42 + $43 + $46 + $47)/4, $18, $52 - $28 -$29 -$30 -$31, $5 - $6, $0}') | tee >(perl -pe 's/ /\t/g' > formike.txt) | ~/mumdex/even_columns ' ' > formike_even.txt
fi

cat denovo_cand.txt | awk '{if ($4 <= '$max_families') print}'  > full_filter_cand.txt
# echo candidates fully filtered $(cat full_filter_cand.txt | wc -l)

# (echo $header ; cat denovo_cand.txt) | ~/mumdex/pivot_table.pl | ~/mumdex/even_columns ' '

show="nFam nInFam nS pNP nSam offset bridge pNB pMed pMax tBC nBC pb sb mcA fcA pcA scA mcB fcB pcB scB mapA mapB supA supB mCA mSA mCB mSB maA faA paA saA maB faB paB saB"

# (echo $header ; cat denovo_cand.txt) | ~/mumdex/pivot_table.pl "nFam nInFam nS pNP" "$show" | ~/mumdex/even_columns ' '


echo
(echo $header ; cat denovo_cand.txt) | ~/mumdex/pivot_table.pl nFam "$show" |
~/mumdex/even_columns ' '

#echo
#(echo $header ; cat denovo_cand.txt) | ~/mumdex/pivot_table.pl nInFam "$show" |
#~/mumdex/even_columns ' '

echo
(echo $header ; cat denovo_cand.txt) | ~/mumdex/pivot_table.pl nS "$show" |
~/mumdex/even_columns ' '

echo
(echo $header ; cat denovo_cand.txt) | ~/mumdex/pivot_table.pl pNP "$show" |
~/mumdex/even_columns ' '
