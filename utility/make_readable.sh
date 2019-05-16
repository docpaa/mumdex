#! /bin/bash

dir="$@"
echo $dir
find $dir -type d | xargs chmod a+rx
find $dir -type f | xargs chmod a+r
