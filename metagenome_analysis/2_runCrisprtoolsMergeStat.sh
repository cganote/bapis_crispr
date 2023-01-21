#!/usr/bin/env bash
#SBATCH -J crasspipeline
#SBATCH -o %x_%j.txt
#SBATCH -e %x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cganote@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=2:00:00
#SBATCH --mem=24G

startSec=$(date +%s) 

## You may need to build crass, libcrispr, and crisprtools on your system.
## I had to compile Crass and Crisprtools configured like so:
# 11-01-22:12:31 ./configure --prefix=/path/where/you/want/crisprtools --with-xerces=/path/to/xerces-c++/3.2.2 --with-libcrispr=/path/to/libcrispr CXXFLAGS=-lstdc++

##############################################################################################################
## Run crass on raw reads mapped against (and filtered for unmaps) Apis mellifera
##############################################################################################################
export LD_LIBRARY_PATH=/path/to/libcrispr/lib:$LD_LIBRARY_PATH
export PATH=/path/to/crass:$PATH

do_the_thing() {
    bam=$1
    basename=${bam##*/}
    rmextra=${basename/unaligned_/}
    sample=${rmextra%.bam}
    left=${sample}_1.fastq
    right=${sample}_2.fastq
    single=${sample}_single.fastq
    ok="new_${sample}.ok"
    if test -f $ok; then
        echo "Already ran $1, no need to redo"
    else
        echo "going to run $1 -> $ok"
	## If you want to retain singletons, send -s to $single instead of /dev/null
        samtools fastq -1 $left -2 $right -0 /dev/null -s /dev/null $bam
        crass ${left} -n 2 -f 2  -o n2_new_${sample}_01 && touch $ok
        crass ${right} -n 2 -f 2  -o n2_new_${sample}_02 && touch $ok
	#crass ${single} -n 2 -f 2  -o n2_new_${sample}_single && touch $ok	
    fi
}

export -f do_the_thing
# I used gnu parallel but you don't have to; just pull do_the_thing out of a function and put it inline.
parallel -j 23 do_the_thing {1} ::: exampledata/unaligned_*.bam

##############################################################################################################
## Run crisprtools on results and extract spacers that fit criteria
##############################################################################################################

crisps=$(find exampledata -name "crass.crispr" | grep "raw")
./crisprtools merge -s $crisps
./crisprtools extract -d crisprtools_merged.crispr > alldrs.fasta
# CRISPRClassify doesn't like fastas, so we'll just pull out the repeat into a txt file. Extension has to be txt =\
grep -v "^>" alldrs.fasta > alldrs.txt
# Run classify on these direct repeats to make sure they are legit
module load r &&  Rscript -e "library(CRISPRclassify); CRISPRclassify::classifyRepeats(\"alldrs.txt\")"
# Clean up the quote marks out of the output
tr -d '"' > alldrs.txt.crclass.clean <alldrs.txt.crclass
for crisp in $crisps; do
    echo $crisp
    #example: n2_new_raw_E4_02/crass.crispr
    nocrass=${crisp%/crass.crispr}
    sample=${nocrass##*n2_new_raw_}
    echo "$nocrass $sample"
    ./crisprtools extract -s -d $crisp | sed -e "s/^>/>${sample}_/g" >> allsamples.fasta
done
perl 3_wrangleDups.pl -r alldrs.txt.crclass.clean -d allsamples.fasta -o allspacersoutput.fasta # default minimum 16

##############################################################################################################
## In our case, crass misses some of the CRISPRs we really care about; pull these out manually and add those
##    spacers to the f file
##############################################################################################################

repeats="CTGTTCCCCGCCCATGCGGGGATGAACCG CGGTTCATCCCCGCATGGGCGGGGAACAG CTGTTCCCCGCCTACGCGGGGATGAACCG CGGTTCATCCCCGCGTAGGCGGGGAACAG GCCGTGGTTTCCCTACCGATTTGTCTATGGTAGCCT AGGCTACCATAGACAAATCGGTAGGGAAACCACGGC"
id=1201 # ideally I would derive this from the end of the allspacersoutput.fasta but for now it's hardcoded
for repeat in $repeats; do
	#grep -h $repeat raw_*.fastq | perl scrapeBombellaRepeats.pl -r $repeat -m 15 > ${repeat}.out
	## Even though this gets me close, and formats things nicely, I still need to do some cleanup manually for these files
	echo  ${repeat}.clean.out
	## The results of manual curation are good but I lost my provenance during that grep step; cleaning up by re-associating which spacer came from where
	while read pattern; do
	    echo "pattern is $pattern"
	    header=$(grep -H -F $pattern raw_*.fastq | perl -ne 'if ( m/raw(.+)\.fastq:/){ $hashy{$1}++;} END { print join("_",  keys %hashy); }')  #>> ${repeat}.clean.out.fast
	    let "id++"
	    #>ID0_E5_02 with repeats  CGGTTCATCCCCGCGTAGGCGGGGAACAG I-E 1.000
	    echo ">ID${id}${header} with repeats  $repeat bapis 1.00" >> allspacersoutput.fasta
	    echo $pattern >> allspacersoutput.fasta
	    let "id++"
	    echo ">ID${id}${header} with repeats  $repeat rev_bapis 1.00" >> allspacersoutput.fasta	    
	    perl revcomp.pl $pattern >> allspacersoutput.fasta
	    echo >> allspacersoutput.fasta	    
	done <${repeat}.clean.out
    done

##############################################################################################################
## Map filtered spacers against bombella genomes and published viromes
##############################################################################################################

cat exampledata/bonillo_isolates.fasta >> virome.fasta
cat exampledata/JAAOBB01.1.fsa_nt >> virome.fasta
genomes=$(find exampledata/ref_genomes -name "Bombella*.fasta")
for genome in $genomes; do
    nopath=${genome##*/}
    noext=${nopath%_GC*.fasta}
    echo "file $genome with tag ${noext}_"
    sed -r -e "s/^>/>${noext}_/g" <$genome >> virome.fasta
done
module load bowtie2/intel/2.4.2
bowtie2-build virome.fasta viromesetc
bowtie2 -x viromesetc --no-unal -a -p 4 -f allspacersoutput.fasta > allspacersoutput_virome_mapped.sam && touch allspacersviromeoutput.ok
grep -v "^@" allspacersoutput_virome_mapped.sam > allspacersoutput_virome_mapped_noheaders.sam


##############################################################################################################
## How do we make sense of all this data? 
##############################################################################################################

    
didItWork=$?
endSec=$(date +%s)

## Don't change lines below:
stop=$(($LINENO - 3))
startTime=$(date -d "@$startSec" +"%m/%d/%y %H:%M:%S")
endTime=$(date -d "@$endSec" +"%m/%d/%y %H:%M:%S")
endFilename=$(date -d "@$endSec" +"%m%d%y_%H_%M_%S")
elapsed=$(($endSec - $startSec))
jobname=$(head -n 2 ${0} | tail -n 1 | perl -pe 's/#SBATCH *-J *(.*)$/$1/') 
me="${jobname:=runlog}"
nospace=${me//[ 	]/}
runlog="${nospace}.${endFilename}"
echo "Ran $me" >> $runlog
echo "Job took" $((elapsed/86400))" days," $(date -d "@$elapsed" -u +"%H hours, %M minutes, and  %S seconds")  >> $runlog 
echo "Starting at $startTime and ending at $endTime" >> $runlog 
echo "Ran on machine $(hostname)" >> $runlog 
echo "Exited with code $didItWork" >> $runlog 
echo "Slurm job ID is ${SLURM_JOB_ID:=NA} " >> $runlog 
echo >> $runlog 
echo "Here's a shameless dump of my contents:" >> $runlog 
head -n $stop ${0} >> $runlog 
