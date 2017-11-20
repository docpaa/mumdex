#! /bin/bash

while [ "$#" -gt 0 ] ; do
   file=$1
   ~/mumdex/bridges2txt ~/analysis/mums/hg19/chrAll.fa $file
   shift
done | perl -ne '
BEGIN {
our %counts;
}
my @fields = split;
if ($fields[0] ne $fields[3]) {
  $type = "t";
} else {
  $type = $fields[6];
}
$invariant = $fields[7];
$absinv = abs($invariant);
if ($absinv < 100) {
  $sinv = $invariant;
} else {
  $sinv = sprintf("%.2g", $invariant);
}

#print $_;
$counts{$type . " " . $sinv} += $fields[13];

END {
  print "stats $#counts\n";
  for my $key (sort {$counts{$a} <=> $counts{$b}} keys %counts) {
    print "$key $counts{$key}\n";
  }
}
'
