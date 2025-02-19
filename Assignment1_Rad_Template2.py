import argparse # Take in arguments
import pyfaidx # Read fasta files
import pdb # Debug
import vcf # Read in vcf files --conda install conda-forge::pyvcf
from Bio.Restriction.Restriction_Dictionary import rest_dict # Library of restriction enzymes so users can enter enzyme names instead of motifs

# Usage: python Assignment1_Rad_Template.py Mzebra_GT3_Chr1.fasta YH_MC_samples_Chr1.vcf SingleRad EcoRI 
# Usage: python Assignment1_Rad_Template2.py Mzebra_GT3_Chr1.fasta YH_MC_samples_Chr1.vcf ddRad EcoRI --RE2 HindIII



"""
Setup: Create a conda environment for this homework and install pyfaidx, vcf, and biopython. Confirm you can import 
each of these modules into python

General Instructions: Below you'll find a code outline for the assignment. Your job is to fill in the missing code
wherever there is a comment starting with "TODO". All of the missing code is contained within functions, which helps
keep individual tasks manageable. Each function already has type-hinted arguments and return variables defined (look up
"type hints in python" if this is confusing), as well as a doc-string with additional instructions. Often, you will
encounter instances where a certain variable is set to "None". These are place-holders, and you will usually want to
replace the "None" with a meaningful definition. 

At the bottom of this script are two sections of code, both commented
out: one called "testing code", the other called "main code". The testing code is meant to be uncommented incrementally,
and is designed to help you check your work as you go. Once you finish filling in all of the missing code segments,
you can uncomment the "main code" section to bring everything together into a working script. 

"""

def my_parse_args() -> argparse.Namespace:
    """
    Task 1:
    Implement a function that uses argparse to parse the command-line options. This function has no arguments itself,
    but should return an argparse.Namespace object named "parsed_args" containing the command-line arguments.
    Set up your parser to take in the following positional arguments (in this order):
        GenomeFile: the path to the .fasta file containing the genome
        VCFFile: the path to the .vcf file with the DNA polymorphisms
        Mode: whether to run SingleRad or ddRad
        RE1: the name (not the motif) of the restriction enzyme
    As well as the following flagged argument:
        RE2: the name (not the motif) of the second restriction enzyme (only used in ddRad mode)

    :return parsed_args: the parsed command-line arguments
    """

    parser = argparse.ArgumentParser(
        description='this script performs either SingleRad or ddRad sequencing on DNA from an input FASTA file and compared'
                    'the results to variants from an input VCF file')
    # TODO: Complete this function as described in the docstring
    parser.add_argument('GenomeFile', type = str, help = 'Path to the .fasta file containing the genome') 
    parser.add_argument('VCFFile', type = str, help = 'Path to the .vcf file with the DNA polymorphisms') 
    parser.add_argument('Mode', type = str, choices = ['SingleRad', 'ddRad'], help = 'What type of RAD to run') 
    parser.add_argument('RE1', type=str, help='Name of restriction enzyme')
    parser.add_argument('-re2','--RE2', type=str, help='Optional argment of name of second restriction enzyme to use in ddRad mode') # optional argument; needs -- or -


    parsed_args = parser.parse_args()
    return parsed_args


def read_fasta(path_to_fasta: str, chromosome: str) -> str:
    """
    Task 2:
    Implement a function to read the fasta file using pyfaidx. This function should return the dna sequence (as a str)
    associated with the specified chromosome.

    :param path_to_fasta: path to the fasta file, as a str
    :param chromosome: which chromosome to read. The default value of 'NC_036780.1' refers to the first chromosome
    :return dna: the dna associated with the specified chromosome, as a str
    """

    # Load the fasta file using pyfaidx
    fasta_file = pyfaidx.Fasta(path_to_fasta)

    dna = str(fasta_file[chromosome][:])  # Convert the sequence to string

    return dna


def find_motifs(dna: str, motif: str): # -> list[int]
    """
    Task 3:
    Implement the find_motifs function to take DNA sequence (string, as returned by your read_fasta function)
    and a Restriction Enzyme motif (string) and will return a list of locations (integers) where the motif is found

    :param dna: DNA sequence to be analyzed
    :param motif: restriction enzyme motif
    :return positions: list of ints of the motif positions
    """
    positions = []
    start = 0 
    while True:
        # pdb.set_trace()
        start = dna.find(motif,start) # .find() gets the first occurence of the motif; start -1 when there is no match
        if start == -1:
            break
        positions.append(start)
        start += len(motif)
    return positions


def run_single_rad(dna: str, re1: str): # -> list[tuple[int, int]]
    """
    Task 4:
    implement a function for performing SingleRad. This function takes in the dna (str) and the name of the restriction
    enzyme (str), and returns a list of tuples where each tuple contains the start index (int) and stop index (int) of a
    site sequenced by SingleRad

    :param dna: DNA sequence to be analyzed
    :param re1: restriction enzyme name
    :return sequenced_sites: a list of tuples, where tuple contains the start and stop index of a sequenced site
    """

    seq_length = 100 # This variable stores how long the sequencing reads are


    # Ensure the enzyme exists in the restriction dictionary
    assert re1 in rest_dict, f"Error: Restriction enzyme {re1} not found in rest_dict"

    # Get the recognition motif for the given enzyme
    motif = rest_dict[re1]['site']

    # Find the cut sites using the find_motifs function
    cut_sites = find_motifs(dna, motif)

    # Generate the sequenced sites with correct range
    sequenced_sites = [(site - seq_length, site + seq_length) for site in cut_sites]

    return sequenced_sites


def run_ddrad(dna: str, re1: str, re2: str): # -> list[tuple[int, int]]:
    """
    Task 5:
    Implement a function for performing ddRad. This function takes in the DNA sequence (str) and the names of the restriction
    enzymes (str and str), and returns a list of tuples where each tuple contains the start index (int) and stop index (int)
    of a site sequenced by ddRad

    :param dna: DNA sequence to be analyzed
    :param re1: first restriction enzyme name
    :param re2: second restriction enzyme name
    :return sequenced_sites: a list of tuples, where each tuple contains the start and stop index of a sequenced site
    """

    # problem with my code is im iterating through each restriction site pair and taking the start and end of it, 
    # but in test code approach is ikely oding the start and stop site of the first of pair and same for the second of the pair


    min_size = 300  # This variable stores the shortest distance for the DNA piece
    max_size = 700  # This variable stores the longest distance for the DNA piece
    seq_length = 100  # This variable stores how long the sequencing reads are

    # TODO: Use an assert statement to verify that RE1 and RE2 are both in the rest_dict provided by biopython
    assert re1 in rest_dict, f"Error: Restriction enzyme {re1} not found in rest_dict"
    assert re2 in rest_dict, f"Error: Restriction enzyme {re2} not found in rest_dict"

    # TODO: Use the rest_dict to look up the motifs associated with the two enzymes
    motif1 = rest_dict[re1]['site']
    motif2 = rest_dict[re2]['site']

    # TODO: Use your find_motifs function to find the sequencing sites for RE1 and RE2
    re1_sites = find_motifs(dna, motif1)
    re2_sites = find_motifs(dna, motif2)

    # Initialize an empty list to store sequenced sites
    sequenced_sites = []

    # Loop through each RE1 site
    for r1 in re1_sites:
        # TODO: Look left of the current site
        left_hits_r1 = [x for x in re1_sites if x < r1]
        left_hits_r2 = [x for x in re2_sites if x < r1]

        # TODO: Check conditions for left side
        if left_hits_r2 and min_size < (r1 - left_hits_r2[-1]) < max_size:
            if not left_hits_r1 or left_hits_r2[-1] > left_hits_r1[-1]:
                # sequenced_sites.append((left_hits_r2[-1], r1))
                sequenced_sites.append((left_hits_r2[-1], left_hits_r2[-1]+seq_length))
                sequenced_sites.append((r1-seq_length, r1))

                # (x, x+100 position) as the first sequenced site
                # and (X - 100, X) as the next sequenced site?

                '''
                In DDRAD sequencing, two different restriction enzymes cut the DNA at sites (x, X).
                The selected fragments lie between these sites and are size-selected.
                Sequencing starts at each cut site and extends outward for a fixed length (e.g., 100 bp).
                This results in reads like (x, x+100) and (X-100, X), rather than sequencing the full x to X region.

                
                '''

        # TODO: Look right of the current site
        right_hits_r1 = [x for x in re1_sites if x > r1]
        right_hits_r2 = [x for x in re2_sites if x > r1]

        # TODO: Check conditions for right side
        if right_hits_r2 and min_size < (right_hits_r2[0] - r1) < max_size:
            if not right_hits_r1 or right_hits_r2[0] < right_hits_r1[0]:
                # sequenced_sites.append((r1, right_hits_r2[0]))
                sequenced_sites.append((r1, r1+seq_length))
                sequenced_sites.append((right_hits_r2[0]-seq_length, right_hits_r2[0]))

    return sequenced_sites




def find_variable_sites(vcf_file_path: str, sequenced_sites):
    """
    Task 6:
    Implement this function to take in a path to a VCF file and a list of sequenced sites (as returned by run_single_rad
    or run_ddrad) and return a list of sites that are both sequenced and contain variation.
    """

    #TODO: read in the VCF file using the VCFReader class
    vcf_obj = vcf.Reader(filename=vcf_file_path)  # Open the VCF file

    sequenced_sites_variable = []  # This list will keep track of all of the sequenced DNA pieces that also contain variation
    
    #TODO: complete the below loop
    for record in vcf_obj:  # Create loop to read through every record
        if record.CHROM != 'NC_036780.1':  # These two lines will make sure that you are only looking at records from Chromosome 1
            continue  # Skip non-chromosome 1 records
        
        # TODO:add an if statement here to ensure that both samples have called genotypes (num_called attribute is your friend)
        if record.num_called != 2:
            continue  # Skip records where both genotypes aren't called

        # TODO:add another if statement to ensure that one genotype is 0/0 and the other is 1/1 (use gt_nums)
        # Ensure one genotype is 0/0 and the other is 1/1 (check the genotype fields)
        genotypes = [sample['GT'] for sample in record.samples]
        if '0/0' in genotypes and '1/1' in genotypes:
            # TODO: if both if conditions pass, use list comprehension to check if the variant's position falls within any sequenced site
            if any(start <= record.POS <= stop for start, stop in sequenced_sites):
                # TODO: if the variant falls within sequenced DNA, add the variant position and chromosome to the list
                sequenced_sites_variable.append((record.POS, record.CHROM))

    return sequenced_sites_variable
# ###TESTING CODE###
# This section of code is provided to help you test your functions as you develop them. Feel free to add additional
# debugging code to this section, but ensure that all code required to run the script remains in the function
# definitions above. As you progress, leave the completed sections uncommented, as some variables are reused. Note
# also that, in much of this testing code, your command-line input (other than file locations) is ignored, but you
# will likely still need to provide a value for each argument to get your code to run.

# Task 1: Parsing Command-line Arguments. Uncomment this section when you have completed the my_parse_args function
# and ensure the args that print match your expectations. This section also checks that the RE1 argument you provide
# is a valid restriction enzyme name from the rest_dict
print('testing my_parse_args function')
args = my_parse_args()
print(f'the following arguments were received: {args}')
if args.RE1 not in rest_dict:
    print(f'no restriction enzyme named {args.RE1} found in rest_dict. Run rest_dict.keys() to see valid names')
print('\n')
#
#
# Task2: Reading the Fasta File. Uncomment this section when you have completed the read_fasta function.
print('testing read_fasta function')
target_chromosome = 'NC_036780.1'
dna_string = read_fasta(path_to_fasta=args.GenomeFile, chromosome=target_chromosome)
print(f'fasta file read. Length of {target_chromosome} is {len(dna_string)}')
print(f'length matches expected length of 45605991: {len(dna_string)==45605991}')
print('\n')
#
# Task 3: Finding Motifs. Uncomment this section when you have completed the find_motifs function.
print(f'testing find_motifs function')
motif_locations = find_motifs('GGTTAAAGATCGGCGAGCCAATGGATCGACGATCA','GATC')
print(f'find_motifs returned {motif_locations}')
print(f'motif locations match expected locations of [7,23,30]: {motif_locations==[7,23,30]}')
print('\n')

# Task 4: Implementing SingleRad. Uncomment this section when you have completed the run_single_rad function. Note that
# this section hard-codes the RE1 argument for consistency, rather than using args.RE1 to get a value from the command
# line
print('testing run_single_rad function')
srad_seq_sites = run_single_rad(dna=dna_string, re1='AanI')
print(f'located {len(srad_seq_sites)} singlerad sites for re1=AanI')
print(f'number of sites matches expected number (19280 sites): {len(srad_seq_sites) == 19280}')
print(f'first sequence site determined to be {srad_seq_sites[0]}')
print(f'first sequence site matches expected site (44450, 44650): {srad_seq_sites[0] == (44450, 44650)}')
print('\n')


# Task 5: Implementing ddRad. Uncomment this section when you have completed the run_ddrad function. Note again that
# the re1 and re2 arguments are hard-coded here for consistency. This section may  take longer to run than others
# due to the complexity of ddrad, but should still run in well under a minute. If it doesn't, check your code
# for infinite loops
print('testing run_ddrad function')
ddrad_seq_sites = run_ddrad(dna_string, re1='AanI', re2='MroI')
print(f'located {len(ddrad_seq_sites)} ddRad sites for re1=AanI, re2=MroI')
print(f'number of sites matches expected number (848 sites): {len(ddrad_seq_sites) == 848}')
print(f'first sequence site determined to be {ddrad_seq_sites[0]}')
print(f'first sequence site matches expected site (91143, 91243): {ddrad_seq_sites[0] == (91143, 91243)}')
print("first seq site:", ddrad_seq_sites[0])
print("425th seq site:", ddrad_seq_sites[424])
print('\n')

# Task 6: Finding Sites with Variation. Uncomment this section when you have completed the find_variable_sites function.
# Again, this section will be slightly slower to run, but should still run in less than a minute
print('testing find_variable_sites function')
variable_srad_seq_sites = find_variable_sites(vcf_file_path=args.VCFFile, sequenced_sites=srad_seq_sites)
print(f'found {len(variable_srad_seq_sites)} singleRad sites with variation')
print(f'number of sites found matches expectation (2483 sites): {len(variable_srad_seq_sites)==2483}')
print(f'{len(100*variable_srad_seq_sites)/len(srad_seq_sites):.3f}% of sites sequenced with single rad had variation')
variable_ddrad_seq_sites = find_variable_sites(vcf_file_path=args.VCFFile, sequenced_sites=ddrad_seq_sites)
print(f'found {len(variable_ddrad_seq_sites)} ddRad sites with variation')
print(f'number of sites found matches expectation (48 sites): {len(variable_ddrad_seq_sites)==48}')
print(f'{len(100*variable_ddrad_seq_sites)/len(ddrad_seq_sites):.3f}% of sites sequenced with ddRad had variation')


### MAIN CODE ###
# Once you think your code is complete, comment out the above testing section, uncomment the below section,
# and run the script.

# args = my_parse_args()
# target_chromosome = 'NC_036780.1'
# dna_string = read_fasta(path_to_fasta=args.GenomeFile, chromosome=target_chromosome)

# if args.Mode == 'SingleRad':
#     print(f'running SingleRad for chromosome {target_chromosome} and restriction enzyme {args.RE1}')
#     seq_sites = run_single_rad(dna_string, args.RE1)
# elif args.Mode == 'ddRad':
#     print(f'running ddRad for chromosome {target_chromosome} and restriction enzymes {args.RE1} and {args.RE2}')
#     seq_sites = run_ddrad(dna_string, args.RE1, args.RE2)

# input_len = len(dna_string)
# n_sites = len(seq_sites)
# sequencing_length = seq_sites[0][1] - seq_sites[0][0]
# sequenced_fraction = sequencing_length*n_sites/input_len
# print(f'{args.Mode} sequencing complete.')
# print(f'\tlength of input sequence: {input_len}')
# print(f'\tnumber of sequencing sites located: {n_sites}')
# print(f'\tpercentage of nucleotides sequenced: {100*sequenced_fraction:.3f}%')

# print('checking for variation within sequenced sites')
# variable_sites = find_variable_sites(args.VCFFile, sequenced_sites=seq_sites)
# n_var_sites = len(variable_sites)
# var_fraction = n_var_sites / n_sites
# print(f'\tnumber of sites with variation: {n_var_sites}')
# print(f'\tfraction of sites with variation: {var_fraction}')

# print('\nAnalysis Complete')


