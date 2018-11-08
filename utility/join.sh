#! /bin/bash

perl -pe '
  print " " if $. > 1;
  chomp;
' 
