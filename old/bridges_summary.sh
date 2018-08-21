#! /bin/bash

exit 0

ref="$1"
shift

while [ "$#" -gt 0 ] ; do
   file=$1
   bridges2txt $ref $file
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
