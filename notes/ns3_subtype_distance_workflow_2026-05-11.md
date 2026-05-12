Create a combined NS3 subtype assignment workbook from the existing NS3 genotype workbook.

Inputs:
- `NS3_Alignments_combined.xlsx`
- `HCV_Subtype_Refs_By_Genome_NA.json`
- the study FASTA directory used for the NS3 genotype workflow

Requirements:
- ignore the excluded quasispecies `RefID`s
- for each retained sequence, only compare against subtype references whose genotype matches the sequence's `BestGT`
- use nucleotide alignment
- skip any alignment shorter than 200 nt

Output rows:
- `RefID`
- `RefName`
- `AccessionID`
- `ClosestGT`
- `ClosestSubtype`
- `ClosestSubtypeDistance`
- `NextClosestSubtypeDistance`
- `AlignedNT`
- `NextClosestSubtypeAlignedNT`

Behavior:
- if multiple reference genomes exist for the same subtype, keep the best hit for that subtype
- determine the closest subtype as the subtype with the lowest uncorrected nucleotide distance
- determine the next closest subtype as the next-best distinct subtype within the same genotype
