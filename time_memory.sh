#!/bin/bash

# Input and output files
input_file="input.fastq"
output_file_trimmomatic="output_trimmomatic.fastq"
output_file_cutadapt="output_cutadapt.fastq"
output_file_fastp="output_fastp.fastq"

# Run Trimmomatic
echo "Running Trimmomatic..."
command_trimmomatic="java -jar trimmomatic-0.39.jar SE $input_file $output_file_trimmomatic TRAILING:3 MINLEN:36"
echo "Running Trimmomatic with the command: $command_trimmomatic"
echo "Trimmomatic Time and Memory Usage:"
/usr/bin/time -v $command_trimmomatic

# Run Cutadapt
echo "Running Cutadapt..."
command_cutadapt="cutadapt -q 20,20 -m 36 -o $output_file_cutadapt $input_file"
echo "Running Cutadapt with the command: $command_cutadapt"
echo "Cutadapt Time and Memory Usage:"
/usr/bin/time -v $command_cutadapt

# Run Fastp
echo "Running Fastp..."
command_fastp="fastp -i $input_file -o $output_file_fastp -q 20 -l 36"
echo "Running Fastp with the command: $command_fastp"
echo "Fastp Time and Memory Usage:"
/usr/bin/time -v $command_fastp

# Output completion message
echo "All processing complete."
