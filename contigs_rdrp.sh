#!/bin/bash

## Rdrp_scan.sh v.1##

# Written by Matheus Cosentino et al 2023.

##############################################################################################################################
# Set functions
help(){

echo "
The folllowing script was built by Matheus Cosentino to use the RdRp-scan workflown developed by Charon et al., 2022
The Script is divided in the following parts
    1. Obtain contigs >600 nts
    2. Blastx with nrdb rdrp_scan (diamond)
    3. Retain all sequences without hits in nr/RdRp_cd
        3.1. Translate ORFs (cut-off 200 aa) using GetORF
        3.2. Remove redundant sequences using CD-Hit    
        3.3. Comparison of HMM-RdRp hits against PDB using Phyre2 server > rdrp-like candidates
    4. Detection of A, B and C motifs using the RdRp motif database (geneious)

For the correct function of the pipeline, install the following mamba enviorments:
    1. Spades (ok)
    2. fastp (ok)
    3. Diamond (ok)
    4. CD-Hit (ok)
    5. Get-ORF (ok)
    6. HMMER (ok)
 
 The packages can be found in mamba enviorment < rdrp_scan >

 Commands:
 bash rdrp_scan.sh --input <dir with libs contigs > --output <dir to be created and add results>
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

echo "#######################################"
echo "## Step 1. Create Output Directories ##"
echo "#######################################"

# Begining of analysis
mkdir -p $output/diamond/ $output/hmmer/

# For each sample
for filename in "$input"/*
do
	# Skip $output/
	if [ $filename == "$output/" ]
	then
		continue
	fi
	# Enter sample diretory
	cd $input

	# Print sample directory name
	echo "Working on -> $filename"
    name=$(echo $filename | sed 's/.*\///')
    echo "$name"

	awk '/^>/ { if (seqlen >= 600) { print header; print seq } header=$0; seq=""; seqlen=0; next }
	 /^[^>]/ { seq = seq $0; seqlen = seqlen + length($0) } END { if (seqlen > 0) 
	{ print header; print seq } }' $name > $name.600nt_contigs.fasta
	mv $name.600nt_contigs.fasta $output/diamond
	
	echo "###################################################################"
	echo "####################  2. Blastx using Diamond #####################"
	echo "###################################################################"
	
	#diamond blastn	
	cd $output/diamond
	diamond blastx -d ~/RdRp-scan/RdRp-scan_0.90.dmnd -q $output/diamond/$name.600nt_contigs.fasta -o $output/diamond/rdrp_$name.txt --min-orf 600 -e le-5 --very-sensitive | 
	awk '{print $1}' rdrp_$name.txt | sort -u > $name.match_contigs.txt |
	grep -Fvf $name.match_contigs.txt $output/diamond/$name.600nt_contigs.fasta > $name.contigs_no_matches.fasta

	echo $name.contigs_no_matches.fasta

	echo "###############################################################################################"
	echo "####################  3. Translate ORFs (cut-off 200 aa) using GetORF #########################"
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
	echo "####################### 4. Remove redundant sequences using CD-Hit ############################"
	echo "##############################################################################################"
	
	#remove duplicated orfs
	cd-hit -i $name.Orfs_600.fasta -o $name.CD_hit_98.fasta -c 0.98 -d 1
	
	#ajust cd-hit output to further use
	awk '{if(NR==1) {print $1} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $name.CD_hit_98.fasta | awk '/^>/{ gsub(" ", "_", $0); } 1' > $name.edited.CD_hit_98.fasta
	mv $output/diamond/$name.edited.CD_hit_98.fasta $output/hmmer/$name.edited.CD_hit_98.fasta

	echo "###############################################################################################"
	echo "####################################### 5. HMM-RdRp hits ######################################"
	echo "##############################################################################################"

	# use hmmer to idd putative rdrp-scan
	cd $output/hmmer
    hmmscan -E 1e-6 --cpu 8 --domtblout $name.hmmscan.txt /home/lddv/RdRp-scan/Profile_db_and_alignments/RdRp_HMM_profile_CLUSTALO.db.h3m $name.edited.CD_hit_98.fasta
	#hmmscan -E 1e-6 --cpu 8 --tblout $name.hmmscan.txt /home/lddv/RdRp-scan/Profile_db_and_alignments/RdRp_HMM_profile_CLUSTALO.db.h3m $name.edited.CD_hit_98.fasta

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

done

date

echo "#######################"
echo "####### The end #######"
echo "#######################"

exit 

done

#hmmscan -E 1e-6 --cpu 8 --domtblout hmmscan_results.txt /home/lddv/RdRp-scan/Profile_db_and_alignments/RdRp_HMM_profile_CLUSTALO.db.h3m contigs.fasta.edited.CD_hit_98.fasta
#grep -v "^#" hmmscan_results.txt | awk '{print $4}' | sort -u > hmmscan_contigs_hits.txt