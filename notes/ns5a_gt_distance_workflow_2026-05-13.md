Create a spreadsheet containing the following fields: RefID and RefName (from the selected worksheet), GenBank accession ID (from the individual FASTA files), the proportion of mismatches (genetic distance) between each individual FASTA file sequence and each GT reference sequence (GT1-GT8) for NS5A, the HCV genotype defined as the closest GT by uncorrected nucleotide distance for that accession ID, and the extent of overlap between the closest matching sequences, defined as the length of a statistically significant nucleotide alignment between the two aligned sequences.

Many sequences in the FASTA files will not contain NS5A and will not be alignable to the reference sequences. These should be skipped.

Steps:

1. Extract the reference GT nucleotide sequences from HCV_GT_RefSeqs.fasta. The relevant information is in the first text string in the header following HCV. The next character (an integer) is the GT; the remaining characters are the gene. For this script, only NS5A is relevant.
2. Extract each of the FASTA nucleotide files from the study FASTA files identified from the worksheet filter.
3. Align each of the sequences in the FASTA files to the reference FASTA nucleotide files and compute the proportion of mismatches.
4. For each sequence in the FASTA files, determine after alignment:
   - the GenBank accession number in the FASTA files
   - the number of nucleotides that were aligned
   - the GT with the lowest proportion of mismatches
   - the proportion of mismatches for each GT
   - the number of nucleotides in the alignment

Create one Excel file for each study FASTA file and place them in a directory called NS5A_Alignments.xlsx so progress can be inspected before later combining them.
