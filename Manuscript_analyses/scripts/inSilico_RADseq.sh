#!/bin/bash

# Load modules
module load bwa/0.7.17
module load samtools/1.9
module load bedtools2/2.24.0

###############################################                                 
echo "---- in silico digest -----"                                 
############################################### 
# Remove blank space between chromosomes
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' \
Macaca_mulatta.Mmul_10.dna.toplevel.fa > mmul_fa_noBreaks.fa

# Remove chr headings
grep -v 'dna:primary_assembly' mmul_fa_noBreaks.fa > mmul_noChrLines.txt

# Remove chr line breaks 
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' mmul_noChrLines.txt > mmul_1Line.txt

# Replace line breaks
sed -e s/GCATG/GCATG'\n'C/g mmul_1Line.txt > mmul_selected_digest.txt

# Filter for inserts that are 250-350
# turn into fasta file
awk '{ if (length($0) < 350) print }' mmul_selected_digest.txt | awk '{ if (length($0) > 250) print }' > mmul_selected_digest_filt.txt

# Get the first 100bp of each read and the last
sed 's/^.*\(.\{100\}\)/\1/' mmul_selected_digest_filt.txt > temp.txt
pr -n:7 -t -T temp.txt | sed 's/^/>/' | tr ":" "\n" | tr -d " " > SphI_digest_100bp_R.txt
cat mmul_selected_digest_filt.txt | colrm 101 > temp.txt
pr -n:7 -t -T temp.txt | sed 's/^/>/' | tr ":" "\n" | tr -d " " > SphI_digest_100bp.txt
rm temp.txt mmul_selected_digest_filt.txt


###############################################
echo "---- mapping and conversion to BED -----"   
###############################################
# Map 
R1_path_input=/inSilico_RADseq/SphI_digest_100bp_R.txt
R2_path_input=/inSilico_RADseq/SphI_digest_100bp.txt
path_genome=Macaca_mulatta.Mmul_10.dna.toplevel.fa
path_out=/inSilico_RADseq/SphI_digest_100bp_PE.sam

# bwa 
bwa mem -t 8 $path_genome $R1_path_input $R2_path_input > $path_out

# sam to bam format 
samtools view -S -b SphI_digest_100bp_PE.sam > SphI_digest_100bp_PE.bam

# bam to bed format 
bedtools bamtobed -i SphI_digest_100bp_PE.bam > SphI_digest_100bp_PE.bed

echo "----- complete -----" 

