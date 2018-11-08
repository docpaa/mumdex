#! /bin/bash

exit 0

if [ "$outdir" = "" ] ; then
    echo no outdir passed | tee -a status.txt | tee >(cat 1>&2)
    exit 1
fi

if [ "$bam" = "" ] ; then
    echo no bam passed | tee -a status.txt | tee >(cat 1>&2)
    exit 1
fi

start=$(date +%s)

cd $outdir

echo clean | tee -a status.txt | tee >(cat 1>&2)
(cd ~/code/mumdex ; make clean > /dev/null)

echo rsync | tee -a status.txt | tee >(cat 1>&2)
rsync -a ~/code/mumdex/ mumdex-code

echo compile | tee -a status.txt | tee >(cat 1>&2)
(cd mumdex-code ; make clean ; make -j mummer merge_mumdex bridges namepair > /dev/null)

node=$(hostname -s)
echo running MUMdex on $node | tee -a status.txt | tee >(cat 1>&2)

export PATH=$outdir/mumdex-code/:$PATH

set -o pipefail

compile=$(date +%s)


if true ; then
echo
echo running mummer | tee -a status.txt | tee >(cat 1>&2)
if mummer -verbose -rcref -l 20 -qthreads 16 -samin ~/analysis/mums/g1k/human_g1k_v37.fasta <(namepair <(samtools view $bam ) minimal) ; then
    echo mummer has succeeded | tee -a status.txt | tee >(cat 1>&2)
else
    echo mummer has failed | tee -a status.txt | tee >(cat 1>&2)
    exit 1
fi

mummer=$(date +%s)

fi

echo
echo running merge_mumdex | tee -a status.txt | tee >(cat 1>&2)
if ! merge_mumdex mumdex 32 16 ; then
    echo merge_mumdex has failed | tee -a status.txt | tee >(cat 1>&2)
    exit 1
else
    echo merge_mumex has succeeded | tee -a status.txt | tee >(cat 1>&2)
fi

merge=$(date +%s)

echo
echo running bridges | tee -a status.txt | tee >(cat 1>&2)
mkdir -p bridges_snp
if (cd bridges_snp ; bridges ../mumdex 16) ; then
    echo bridges succeeded | tee -a status.txt | tee >(cat 1>&2)
else
    echo bridges has failed | tee -a status.txt | tee >(cat 1>&2)
    exit 1
fi

bridges=$(date +%s)

echo timing compile $((compile - start)) mummer $((mummer - compile)) merge $((merge - mummer)) bridges $((bridges - merge)) all $((bridges - start))


echo
(cd mumdex-code ; make clean > /dev/null)
echo all done | tee -a status.txt | tee >(cat 1>&2)
exit 0
