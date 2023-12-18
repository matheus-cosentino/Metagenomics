#!/bin/bash

## Rdrp_scan.sh v.1##

# Written by Matheus Cosentino et al 2023.

##############################################################################################################################
# Set functions
help(){

echo "

The folllowing script was built by Matheus Cosentino to use the RdRp-scan workflown developed by Charon et al., 2022
The Script is divided in the following parts
    0. Commands verification (ok)
    1. Size and quality trimming of Raw Reads (ok)
    2. De novo Contig Assembly and Size filtering (ok)
    3. Blastx using DIamond and retain all sequences with and without hits in nr/RdRp_cd
    5. Translate ORFs (cut-off 200 aa) using GetORF
    6. Remove redundant sequences using CD-Hit
    7. Comparison of HMM-RdRp hits against PDB using Phyre2 server

	For the correct function of the pipeline, install the following mamba enviorments:
    1. Spades
    2. fastp
    3. Diamond
    4. CD-Hit
    5. Get-ORF
	6. HMMER
	7. trimmal
 
 The packages can be found in mamba enviorment < rdrp_scan >

 Commands:
 bash rdrp_scan.sh --input <dir with libs> --output <dir to be created and add results> --rdrp <dir with rdrp metadata> --threads <number of threads> --adapters <file with fasta adapters>
	bash ~/RdRp-scan/rdrp_scan.sh --input /home/lddv/ssd2/Matheus/SCRIPT_RDRP_SCAN/backup --output /home/lddv/ssd2/Matheus/Rdrp_test --threads 5 --adapters ~/RdRp-scan/ADAPTERS/adapters.fa
"
}

version(){
	echo "Rdrp Scan v0.0.1"
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

########
# Add variables that will be used
# directory with imput

echo "###################################################"
echo "##Step 1. Size and quality trimming of Raw Reads ##"
echo "###################################################"

# Go to libraries root directory and set $output/*/
cd $input

# Begining of analysis

mkdir -p $output/qc_report/raw/ $output/qc_report/filtered/ $output/diamond/ $output/hmmer/ $output/megahit/ 

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

	#make dir of individual samples
	mkdir $output/qc_report/raw/$name

	#move fastqc files to new directory
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

	#make dir of individual samples

	mkdir $output/qc_report/filtered/$name
	mv *.html $output/qc_report/filtered/$name
	mv *.zip $output/qc_report/filtered/$name
	
	# Concatenate unpaired reads
	cat trim.u.*fastq > trim.uni.$R1
	rm trim.u.*

	# Get the names of qc reads
	U=trim.uni.$R1
	R1=trim.p.$R1
	R2=trim.p.$R2

	echo "##################################################################"
	echo "#######  2. De novo Contig Assembly and Size filtering ###########"
	echo "##################################################################"

	#de novo assembly
	megahit -1 $R1 -2 $R2 -o $output/megahit/$name
	#spades.py --meta -1 $R1 -2 $R2 -s $U -o $output/metaspades -t $threads
	mkdir $output/diamond/$name
	cp $output/megahit/$name/final.contigs.fa $output/diamond/$name/$name.contigs.fasta
	#cp $output/metaspades/contigs.fasta $output/diamond/$name.contigs.fasta

	# Compress fastq	
	echo "###################Compressening Librarys##########################"
	mv trim.* $output/qc_report/filtered/$name
	cd $sample
	gzip *fastq
	echo "############################ Done #################################"
	
	# filtering by size
	cd $output/diamond/$name
	awk '/^>/ { if (seqlen >= 600) { print header; print seq } header=$0; seq=""; seqlen=0; next }
	 /^[^>]/ { seq = seq $0; seqlen = seqlen + length($0) } END { if (seqlen > 0) 
	{ print header; print seq } }' $name.contigs.fasta > $name.600nt_contigs.fasta

	echo "###################################################################"
	echo "####################  3. Blastx using Diamond #####################"
	echo "###################################################################"

	#diamond blastn
	diamond blastx -d ~/RdRp-scan/RdRp-scan_0.90.dmnd -q $name.600nt_contigs.fasta -o rdrp_$name.txt --min-orf 600 -e le-5 --very-sensitive | 
	awk '{print $1}' rdrp_$name.txt | sort -u > $name.match_contigs.txt |
	grep -Fvf $name.match_contigs.txt $name.600nt_contigs.fasta > $name.contigs_no_matches.fasta

	echo "###############################################################################################"
	echo "####################  4. Translate ORFs (cut-off 200 aa) using GetORF #########################"
	echo "##############################################################################################"

	# Define the list of genetic code tables
	genetic_code_tables="1 3 4 5 6 11 16"

		# Loop through each genetic code table and run getorf
		for table in $genetic_code_tables; do
			echo "################## Finding ORFS for "$table" genetic code ##########################"
			output_file="ORFs_Frame${table}_"$name"_min600.fasta"
    		getorf $name.contigs_no_matches.fasta -minsize 600 -table "$table" -find 0 -outseq "$output_file"
    		echo "################### ORFS found for "$table" genetic code ##########################"
		done

	# Concatenating multiple files in only one fasta and removing originals from getorf
    cat ORFs_Frame* > $name.Orfs_600.fasta
    rm ORFs_Frame*

	echo "###############################################################################################"
	echo "####################### 5. Remove redundant sequences using CD-Hit ############################"
	echo "##############################################################################################"
	
	#remove duplicated orfs
	cd-hit -i $name.Orfs_600.fasta -o $name.CD_hit_98.fasta -c 0.98 -d 1
	
	#ajust cd-hit output to further use
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $name.CD_hit_98.fasta | awk '/^>/{ gsub(" ", "_", $0); } 1' > $name.edited.CD_hit_98.fasta
	mkdir $output/hmmer/$name
	cp $output/diamond/$name/$name.edited.CD_hit_98.fasta $output/hmmer/$name/$name.edited.CD_hit_98.fasta

	echo "###############################################################################################"
	echo "####################################### 6. HMM-RdRp hits ######################################"
	echo "##############################################################################################"

	# use hmmer to idd putative rdrp-scan
	cd $output/hmmer/$name
    hmmscan -E 1e-6 --cpu 8 --domtblout $name.hmmscan.txt /home/lddv/ssd2/Matheus/Metagenomics_pipe_2023/RDRP_SCAN/Profile_db_and_alignments/RdRp_HMM_profile_CLUSTALO.db.h3m $name.edited.CD_hit_98.fasta
	
	#get a list with sequences matching
	grep -v "^#" $name.hmmscan.txt | awk '{print $4}' | sort -u > $name.hmmscan_contigs_hits.txt
	
	#subset from original fasta
	grep -A1 -Fwf $name.hmmscan_contigs_hits.txt $name.edited.CD_hit_98.fasta > putativeRdrp_$name.fasta
	awk '/^>/ { if (header[$0]++) skip=1; else skip=0 } !skip' putativeRdrp_$name.fasta > unique_putativeRdrp_$name.fasta
	sed 's/--$//' unique_putativeRdrp_$name.fasta > Final_Rdrp_$name
	
	# deleting temporary files
	rm unique_putativeRdrp_$name.fasta
	rm putativeRdrp_$name.fasta 
	rm $name.hmmscan_contigs_hits.txt 
	#rm $name.edited.CD_hit_98.fasta

	cd $input

done

date

echo "#######################"
echo "####### The end #######"
echo "#######################"

exit 