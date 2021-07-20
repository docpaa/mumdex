#! /bin/bash

# Configuration and code details
home=/gpfs/commons/home/pandrews
code=$home/code/smash
bin_dir=$home/bins
ref=$home/marvel/hg19/chrAll.fa

# Run details
data_dir=$1
samples=$(ls $data_dir | cut -d "_" -f 1 | sort | uniq)
n_samples=$(echo $samples | wc -w)
n_fastqs=$(ls $data_dir/*.gz | wc -l)
n_threads=12

echo Working on $n_samples samples $samples with $n_fastqs fastqs using $n_threads threads

bins=20000,50000,100000,200000,500000,1000000

for sample in $samples ; do
    r1=$(ls $data_dir/$sample*.R1.fastq.gz)
    r2=$(ls $data_dir/$sample*.R2.fastq.gz)
    echo Running smash on sample $sample and fastqs $r1, $r2
    if [ ! -e $sample ] ; then
        mkdir -p $sample
        (
            echo '#! /bin/bash'
            echo $code/src/smash_marvel.sh $code $ref $bin_dir $bins $sample '"'$r1'"' '"'$r2'"' $n_threads '&& echo Success' 
        ) | sbatch --job-name=$sample --mem=40G -N 1 --cpus-per-task=$n_threads -o $sample/${sample}.slurm.out.txt -e $sample/${sample}.slurm.err.txt
    fi
done





 
