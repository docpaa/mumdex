#! /bin/bash

set -e -u -E -o pipefail
trap "exit 1" int

while true ; do
    eval "$@" | less -K 2>&1 || break
    sleep 1 || break
done

