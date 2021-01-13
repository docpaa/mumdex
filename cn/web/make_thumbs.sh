#! /bin/bash

data_dir=$1
dataset=${data_dir%/*}
dataset=${dataset##*/}
echo data_dir is $data_dir dataset is $dataset
main_dir=$PWD

# abspos by chromosome
for genome in m38 hg19 ; do
    cat $genome.fa.fai |
        grep -v -e _ -e M |
        awk 'BEGIN {start=0;} { print $1, start, start + $2; start=start + $2; }
END { print "all", 0, start; }  
' > $genome.starts.txt     
done

# create thumbnails
chrs=$(echo all chr{{1..22},{X,Y}}) &&
# chrs=$(echo all chr1) &&
rm -f $data_dir/samples.txt &&
rm -f pngs.txt &&
set pipefail && 
ls -dr $data_dir/*/*.txt |
while read file ; do
    sample_dir=${file%/*}
    thumbs_dir=$sample_dir/thumbs
    mkdir -p $thumbs_dir
    chmod 711 $thumbs_dir
    genome=$(cat $file | grep chr21 > /dev/null && echo hg19 || echo m38)
    # echo $genome > $sample_dir/genome.txt
    chrpos=$(
        head -n 1 $file |
        cut -f 1,2 |
        perl -pe 's/\t/,/g ; s/bin.chrom/bin.chrom:c/'
    )
    sample=$(
        echo ${file##*/} |
        cut -d . -f 1,2
    )
    tumor=${sample%%.*}
    # sample=${sample#*.}

    echo $file $sample $tumor $chrpos $genome
    
    for scale in log linear atan ; do
        for range in common sample ; do
            for color in red blue ; do
                for chr in $chrs ; do
                    log_opt="--$scale"
                    colors=$([ $color = blue ] && echo " --colors 1,0"])
                    [ $scale = atan ] && [ $range = sample ] && log_opt="--atanS"
                    command=$(
                        echo ~/mumdex/ggraph $log_opt $colors --output $(
                            if [ $range = common ] || [ $chr != all ] ; then
                                echo -e --initial $(
                                    cat $genome.starts.txt |
                                        perl -ne "print if /$chr\s/" |
                                        awk '{print $2,$3}'
                                     ) $(
                                    if [ $range = sample ] ; then
                                        echo X X
                                    else
                                        if [ $scale = "log" ] ; then
                                            low=0.1 &&
                                                high=20 &&
                                                echo $(math "log($low/2)/log(10)") \
                                                     $(math "log($high/2)/log(10)")
                                        elif [ $scale = "atan" ] ; then
                                            echo -0.01 1.01
                                        else
                                            echo 0 4.1
                                        fi
                                    fi
                                     )
                            else
                                echo ""
                            fi
                             ) cn $main_dir/$genome.fa \
                                 $(( [ $dataset = breast ] || [ $dataset = micro ] ) &&
                                       echo $chrpos,lowratio,seg.mean.LOWESS $file ||
                                           echo abspos,ratio,seg_ratio $file)
                           )
                    
                    # echo "   " $chr $scale $range $color

                    id=${sample}_${chr}_${range}_${scale}_${color}
                    thumb_dir=$sample_dir/thumbs/$id
                    mkdir -p $thumb_dir
                    chmod 711 $thumb_dir
                    (cat > $thumb_dir/command.sh <<EOF
cd $thumb_dir
if [ -e $id.png ] ; then exit 0 ; fi
rm -f cn.*.{png,xpm,pdf}
(export DISPLAY=:19 &&
$command 2>&1 &&
mv cn.1.png $id.png) |
grep -v -e ^Saved -e ^Converted
EOF
                    )
                    echo $thumb_dir/$id.png bash $thumb_dir/command.sh >> pngs.txt
                done
                
            done
        done
    done
    echo $tumor $sample >> $data_dir/samples.txt
done

echo samples file
cat $data_dir/../*/samples.txt > $data_dir/../samples.txt


exit 0
echo making pngs

(
    Xvfb :19 &
    sleep 2
    cat pngs.txt | while read png command ; do
        echo $command
    done | xargs -I CMD -P 32 bash -c CMD
    killall Xvfb
)

cat pngs.txt |
    while read png command ; do
    if [ ! -e $png ] ; then
        echo $command
    fi
    done | tee $data_dir/missing.txt
rm -f pngs.txt starts.txt ggraph.cfg.new
[ -s $data_dir/missing.txt ] || rm $data_dir/missing.txt
    
