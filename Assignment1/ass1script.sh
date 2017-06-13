#!/bin/bash
#Bash shell will run this script
#This will perform multiple alignment on each provided fasta file in the input directory

echo "Running Multiple Alignment script"

#Path to directory with fasta files
input=/home/jbergenstrahle/BioMedSchool/Assignment1/yeast_genes/*.fasta

#Make directory if needed (-p)
output="/home/jbergenstrahle/BioMedSchool/Assignment1/output/"
mkdir -p /home/jbergenstrahle/BioMedSchool/Assignment1/output

for file in $input
do
	echo "Processing $file ..."
	basename=`basename $file`
	/home/jbergenstrahle/Downloads/muscle -in $file -out $output$basename".afa"
	echo "$file aligned"
done

