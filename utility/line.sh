#! /bin/bash

line=${1:-1}
shift
cat "${@:-/dev/stdin}" |
head -n $line |
tail -n 1
