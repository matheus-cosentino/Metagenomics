#!/bin/bash

## Rdrp_scan.sh v.1##

# Written by Matheus Cosentino et al 2023.

##############################################################################################################################
# Set functions
help(){

echo "

The folllowing script was built by Matheus Cosentino to analyze Viral Metagenomics data.
The Script is divided in the following parts
    0. Commands verification (ok)
    1. Size and quality trimming of Raw Reads
    2. Mapping to Host Genome
    2. De novo Contig Assembly and Size filtering by Megahit
    3. Diamond Blastx against Viral Database
    
 The packages can be found in mamba enviorment < rdrp_scan >

 Commands:
 bash Cosentino_metagen.sh --input <dir with libs> --output <dir to be created and add results> --rdrp <dir with rdrp metadata> --threads <number of threads> --adapters <file with fasta adapters>
	bash ~/Metagenomics_pipe_2023/Cosentino_metagen.sh --input /home/lddv/RdRp-scan/CpPRHant --output /home/lddv/ssd2/Matheus/Metagen_test --threads 5 --adapters ~/RdRp-scan/ADAPTERS/adapters.fa --reference /home/lddv/Metagenomics_pipe_2023/GCF_022682495.1_HLdesRot8A_genomic.fna
"
}

version(){
	echo "Metagenomivcs v0.0.1"
}

##############################################################################################################################
# Validate arguments

# Verify help function
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	help
	exit
fi

# Verify Version function
if [ -z "$1" ] || [[ $1 == -v ]] || [[ $1 == --version ]]; then
	version
	exit
fi

# Verifying arguments in bash comand line
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        echo $1 $2
   fi
   shift
done

# Add variables that will be used
# directory with imput

echo "##############################"
echo "##Step 0. Veryfing Arguments##"
echo "##############################"

if [[ -d $input ]] 
then 
	echo "Root directory for samples indentified ->  $input"
else 
	echo "Samples directory not found. Please specify."
	exit
fi

# directory for output
if [[ -d $output ]] 
then 
	echo "Output directory indentified ->  $output"
	echo "Please specify another to avoid overwriting previous analysis."
	exit
else 
	echo "Samples directory created -> $output"
	#mkdir $output
fi

# trimming command (cmd info)
if [ -z "$trim" ]
then
      echo "Number of bases to trim from the beggining of read not specified, using the default (30bp)"
	  trim=30
else
      echo "Number of bases to trim from the beggining of read: $trim"
fi

# number of threads to use
if [ -z "$threads" ]
then
      echo "Number of threads not specified, starting with the default (1)"
	  threads=1
else
      echo "Number of threads available for processing: $threads"
fi

# number of minimum read length
if [ -z "$minimum_length" ]
then
      echo "Minimum read length not specified, using the default (50bp)"
	  minimum_length=50
else
      echo "Minimum read length: $minimum_length"
fi

# file path with adapters
if [[ -f "$adapters" ]]
then
      echo "Trimmomatic adapters fasta file identified ->  $adapters"
else
    echo "Trimmomatic adapters fasta file not found. Please specify."
    exit
fi

# file path with host genome
if [[ -f "$reference" ]]
then
      echo "Host genome file identified ->  $reference"
else
    echo "Host genome file not found. Please specify."
    exit
fi

# file path with db
if [[ -d "$database" ]]
then
      echo "Database to viral taxonomy iddentiffied ->  $database"
else
    echo "Database to viral taxonomy not iddentiffied. Please specify."
    exit
fi

########
# Add variables that will be used
# directory with imput

echo "###################################################"
echo "##Step 1. Size and quality trimming of Raw Reads ##"
echo "###################################################"

# Go to libraries root directory and set $output/*/
cd $input

# Begining of analysis

mkdir -p $output/qc_report/raw/ $output/qc_report/filtered/ $output/diamond/ 

# For each sample
for sample in */
do
	# Skip $output/
	if [ $sample == "$output/" ]
	then
		continue
	fi

	# Enter sample diretory
	cd $sample

	# Print sample directory name
	echo "Working on -> $sample"

	# Decompress data
	gunzip *gz
	
	# Check raw data files
	N=$(ls *fastq | wc -l)
	if [ $N == 2 ]
	then
		R1="$(ls *fastq | head -1)"
		R2="$(ls *fastq | tail -1)"
	else
		echo "Please provide only two fastq files (R1 and R2) for each sample."
		echo "Only files with the .fastq extension are used."
		echo "Different number of files found in directory $sample."
		exit
	fi
	
	# Get sample name. It should never include a _ character (as it is used as separator)
	name=$(echo "$R1" | grep -oP "^[A-Za-z0-9]+")
	
	# QC report for raw data 
	fastqc -q -t $threads *fastq
	mv *zip $output/qc_report/raw/
	mv *html $output/qc_report/raw/
		
	# Filter data with fastp
	echo "Performing strict QC..."
	if [ $trim -eq 0 ]
	then
		trimmomatic PE -threads $threads -phred33 $R1 $R2 trim.p.$R1 trim.u.$R1 trim.p.$R2 trim.u.$R2 ILLUMINACLIP:$adapters:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:$minimum_length
	else
		trimmomatic PE -threads $threads -phred33 $R1 $R2 trim.p.$R1 trim.u.$R1 trim.p.$R2 trim.u.$R2 ILLUMINACLIP:$adapters:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 HEADCROP:$trim MINLEN:$minimum_length
	fi

	# QC report for filtered data
	echo "Performing QC report for filtered data"
	fastqc -q -t $threads trim*fastq
	mv *html $output/qc_report/filtered/
	mv *zip $output/qc_report/filtered/
	
	# Concatenate unpaired reads
	cat trim.u.*fastq > trim.uni.$R1
	rm trim.u.*

	# Get the names of qc reads
	Ut=trim.uni.$R1
	R1t=trim.p.$R1
	R2t=trim.p.$R2
    	
	echo "#######################################################################"
	echo "#################  2. Host mapping and filtering  #####################"
	echo "#######################################################################"

    #Filter genomic host DNA
	minimap2 -a -t $threads -x sr $reference $R1t $R2t > $name.sam
	#sort sam file
	samtools sort -@ $threads $name.sam -o sorted_$name.sorted.sam
	#convert to bam
	samtools view -@ $threads -bS sorted_$name.sorted.sam  > $name.sorted.bam
	#remove reads not mapped
	bamToFastq -i $name.sorted.bam -fq unmapped_1_$name.fq -fq2 unmapped_2_$name.fq

	echo "##################################################################"
	echo "#################  3. De novo Contig Assembly  ###################"
	echo "##################################################################"
	megahit -1 unmapped_1_$name.fq -2 unmapped_2_$name.fq -o $output/megahit
	
	# Compress fastq	
	echo "####################### Cleaning files ############################"
	rm *bam
	rm *sam
	rm trim*
	rm unmapped*
	echo "###################Compressening Librarys##########################"
	gzip *fastq
	
	echo "####################### Adjusting Assembly ############################"
	cd $output/megahit
	mv final.contigs.fa $name.contigs.fas
	cp $name.contigs.fas $output/diamond/
	echo "############################ Done #################################"

	echo "###################################################################"
	echo "#################  4. Diamond against ViralDB ####################"
	echo "##################################################################"

	cd $output/diamond/
	#blast against viral database
	#diamond blastx -query  -db $database -out $name.Diamondout.txt -evalue 1e-5 -num_threads $threads -outfmt 6
	blastn -query $name.contigs.fas -db $database/ref_viruses_rep_genomes -out $name.Diamond.txt -evalue 1e-5 -num_threads 10 -outfmt 6

    #blast and krona
	#cd $output/diamond/
	#diamond blastx -d $database -q $name.contigs.fasta -o rdrp_$name.txt --min-orf 600 -e le-5 --very-sensitive

done

date

echo "#######################"
echo "####### The end #######"
echo "#######################"

exit

