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
echo This has been tested on Linux and Mac.
echo

# Get and compile ggraph
if [ -e mumdex/ggraph ] ; then
    echo mumdex directory and ggraph program already exists - skipping 1>&2
else
    wget http://mumdex.com/mumdex.zip
    unzip mumdex.zip
    rm mumdex.zip
    cd mumdex
    make -j 4 ggraph
    if [ $? != 0 ] ; then
        echo Problem compiling ggraph - \
            please try to copile manually before continuing - quitting 1>&2
        exit 1
    fi
    cd ../
fi

# Get sample data
data="{m,f,d,s}.txt"
if [ -e s.txt ] ; then
    echo sample data already exists - skipping 1>&2
else
    eval wget http://mumdex.com/ggraph/data/$data
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
        eval wget http://mumdex.com/ggraph/config/$genome/{$config,$gzip}
        echo unzipping $gzip
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
    echo ./mumdex/ggraph cn $genome.fa abspos,ratio,seg "$data"
    if [ ! -e $genome.fa.bin/ref.seq.bin ] ; then
        echo
        echo Note the first time the command runs it will take a short time \
            to create a binary reference cache file
    fi
done
