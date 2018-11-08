#! /bin/bash

seq=$1

echo $seq |
perl -ne '
chomp;
tr/ACGT/TGCA/;
my @seq = reverse split //, $_;
print join("", @seq) . "\n";
'
