#! /bin/bash

defwidth=2560
defheight=1600
altwidth=defwidth
altheight=defheight
altwidth=1920
altheight=1200
#altwidth=2560
#altheight=1600
banner=25

if [ "$1" = "-p" ] ; then
    preload=1
    shift
fi

family=$1
chr=$2
position=$3
message=$4
width=${5:-$altwidth}
height=${6:-$altheight}

height=$((height - banner))
defheight=$((defheight - banner))
width=$((width - edge))
defwidth=$((defwidth - edge))


low=${position%-*}
high=${position#*-}
if [ "$low" = "$high" ] ; then
    low=$((low - 1))
fi

people="father mother self sibling"
family_dir="/data/safe/paa/analysis/mums-new/output/families"

edge=$((10 * width / defwidth))
xoff=$((width / 2))
yoff=$((banner + height / 2))

xterm_width=$((204 * width / defwidth))
xterm_height=$((55 * height / defheight))

x=$edge
y=$((banner + edge))

pids=""
for ind in $people ; do
    args="'$chr $low $high' $family_dir/$family/$ind/rmdup.bam"
    if [ "$preload" = 1 ] ; then
        eval gav merged $args 2> /dev/null > /dev/null &
    else
        xterm -geometry ${xterm_width}x$xterm_height+$x+$y -T "$family $chr $position $ind $message" -e "gav $args ; sleep 6000" &
    fi
    pids="$pids $!"
    if [ $x -eq $xoff ] ; then
        x=$edge
        y=$yoff
    else
        x=$xoff
    fi
done

xoff=$((width/2-190))
yoff=$((banner+height/2-2*edge))

if [ "$preload" = 1 ] ; then
    wait
else
    sleep 1
    xterm -geometry 30x1+$xoff+$yoff -bg red -T "close to quit" -e 'sleep 6000'
    kill -9 $pids
fi

exit 0



