#! /bin/bash

# Genome selection
genomes=$@
if [ "$genomes" != "hg19" ] &&
    [ "$genomes" != "hg38" ] &&
    [ "$genomes" != "hg19 hg38" ] ; then
    echo Required command line arguments must be exactly \
        hg19, hg38, or hg19 hg38 1>&2
    exit 1
fi

echo This program will download, install, and setup G-Graph for you.
echo It will also download the "genome(s)" for \
    $genomes, annotations, and sample data.
echo Please contact Peter Andrews at paa@drpa.us if you encounter difficulties.
echo This has been tested on Linux, Mac, and Windows under Cygwin.
echo

# Special Cygwin setup on Windows
STARTX=""
if [ "$OSTYPE" = cygwin ] ; then
    STARTX=startxwin
    if [ ! -e /bin/wget.exe ] ; then
        /cygdrive/c/Users/$USER/Downloads/setup-x86_64.exe -q -d -n -P wget
    fi

    while [ ! -e /bin/wget.exe ] ; do
        echo Waiting for wget install to finish
        sleep 2 
    done
    sleep 1

    if [ ! -e /usr/local/bin/apt-cyg ] ; then
        wget https://raw.githubusercontent.com/transcode-open/apt-cyg/master/apt-cyg
        chmod +x apt-cyg
        mv apt-cyg /usr/local/bin
    fi

    apt-cyg install git perl-libwww-perl unzip make gcc-g++ libX11-devel xinit ImageMagick xorg-x11-fonts-Type1 libgsl-devel python

fi

# Get and compile ggraph
if [ -e mumdex/ggraph ] ; then
    echo mumdex directory and ggraph program already exists - skipping 1>&2
else
    if [ -e mumdex/ ] ; then
        echo mumdex directory already exists - skipping 1>&2
    else
        url=http://mumdex.com/mumdex.zip
        echo Downloading $url
        lwp-request -m GET $url > mumdex.zip
        echo Unzipping mumdex.zip
        unzip mumdex.zip > /dev/null
        rm mumdex.zip
    fi
    cd mumdex
    echo Compiling mumdex/ggraph and mumdex/x11plot
    make -j 4 ggraph x11plot > /dev/null
    if [ $? != 0 ] ; then
        (
            echo Problem compiling ggraph and x11plot
            echo Please try to copile manually before continuing
            echo You probably need a more recent compiler "(gcc 4.9 or later)"
            echo Quitting
        ) 1>&2
        exit 1
    fi
    cd ../
fi

# Get sample data
data="{m,f,d,s}.txt"
if [ -e s.txt ] ; then
    echo sample data already exists - skipping 1>&2
else
    for file in $(eval echo $data) ; do
        url=http://mumdex.com/ggraph/data/$file
        echo Downloading $url
        lwp-request -m GET $url > $file
    done
fi

# Get genomes
for genome in $genomes ; do
    gzip=$genome.fa.gz
    fasta=$genome.fa
    bin=$genome.fa.bin
    if [ -e $bin ] ; then
        echo Genome files for $genome already exist - skipping 1>&2
    else
        config="{knownGene,knownIsoforms,kgXref,cytoBand}.txt"
        for file in $(eval echo $config $gzip) ; do
            url=http://mumdex.com/ggraph/config/$genome/$file
            echo Downloading $url
            if [ $file = $gzip ] ; then
                echo Please be patient for this one
            fi
            lwp-request -m GET $url > $file
        done
        echo Unzipping $gzip
        gunzip $gzip
        mkdir $bin
        eval mv $config $bin
    fi
done

echo 
echo Setup all done!
for genome in $genomes ; do
    echo
    echo You can now run G-Graph with $genome using the following command:
    echo
    echo $STARTX ./mumdex/ggraph cn $genome.fa abspos,ratio,seg "$data"
    if [ ! -e $genome.fa.bin/ref.seq.bin ] ; then
        echo
        echo Note the first time the command runs it will take a short time \
            to create a binary reference cache file
    fi
    if [ $genome = hg38 ] ; then
        echo 
        echo Note the sample data is processed for hg19 \
            but it will still display "(slightly incorrectly)" for hg38
    fi
done
