#! /bin/bash

# Process arguments
usage="usage: find_dupes.sh [-t n_threads] dir ..."
n_threads=10
if [ "$1" = -t ] ; then
    n_threads=$2
    shift 2
fi
if [ "$#" -lt 1 ] ; then
    echo Not enough arguments 1>&2
    echo $usage 1>&2
    exit 1
fi

echo recording md5 sums for files in $@ using $n_threads threads

# get md5s and group files
find "$@" -type f -print0 |
    xargs -0 -n 10 -P $n_threads md5sum |
    cat |
    perl -ne '
        BEGIN { our %items }
        chomp;
        my ($md5, $file) = /(\w+)\s+(.+)/;
        push @{$items{$md5}}, $file;
        END {
            for my $md5 (keys %items) {
                print $md5 . "\t" . scalar(@{$items{$md5}}) . "\t" .
                join( ";", @{$items{$md5}}) . "\n";
            }           
        }
' |
    sort -k2n



