#! /bin/bash

ls | while read entry ; do
    echo $(du -ks "$entry") $(stat -c %y "$entry" | cut -d . -f 1)
done |
sort -k1n |
even_columns ' '
