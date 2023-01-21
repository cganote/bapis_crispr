#!/usr/bin/env bash 
#SBATCH -J indivqueenmapapimelraw
#SBATCH -o %x_%j.txt
#SBATCH -e %x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cganote@iu.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=2:00:00
#SBATCH --mem=10G

#Load any modules that your program needs
module unload gcc
module load gcc/9.3.0
module load bowtie2
module load samtools

# Change this to a directory containing raw reads for your metagenome
readsdir=exampledata/01_all_raw_reads

# Build the bee index, point to your bee genome fasta file
bowtie2-build exampledata/Apismellifera_DH4_GCF_003254395.2.fasta apimel

# Use whatever scheme your reads are named to grab the left reads first 
for leftread in $readsdir/*_R1_001.fastq.gz; do
    # Change left to right
    rightread=${leftread/_R1_/_R2_}
    # Strip out the path
    removePath=${leftread#${readsdir}/}
    # For my scheme, I wanted to remove everything from the file name after the first _
    sample=${removePath%%_*}
    
    echo "$leftread $rightread $removeTrim $sample"
    # Skip bowtie if it was already run on this sample name
    if test -f  ${sample}_rawbam.ok; then
        echo "already see a file called ${sample}_rawbam.ok, skipping since it already is run."
    else
	## I think -f12 (AND operation) is a bit more lenient than what I was thinking, -f4 -f8 (OR operation)
	echo "did not find $sample, running bowtie" 
	bowtie2 -x apimel -p 24 -1 $leftread -2 $rightread | samtools view -b -f12 > unaligned_raw_${sample}.bam  && touch ${sample}_rawbam.ok        
    fi


done


