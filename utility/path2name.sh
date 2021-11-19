#! /bin/bash

for file in "$@" ; do
    nodown="${file//..\//}"
    nospaces="${nodown// /.}"
    echo ln -sf "$file" "${nospaces//\//.}"
done

    
