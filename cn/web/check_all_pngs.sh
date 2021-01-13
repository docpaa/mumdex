#! /bin/bash

Xvfb :19 &
( echo echo dummy
  for dir in site/data/*/*/*/thumbs/* ; do
      id=${dir##*/}
      png=$dir/$id.png
      if [ ! -e $png ] ; then
          echo making $png 1>&2
          echo bash $dir/command.sh
      fi
  done
) | xargs -I CMD -P 32 bash -c CMD

killall Xvfb

find site -name 'cn*.pdf' | xargs rm
find site -name 'cn*.xpm' | xargs rm
find site -name 'ggraph.cfg.new' | xargs rm
