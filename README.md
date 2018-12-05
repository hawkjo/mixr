# MIXR: Mismatching Isoform eXon Remover

Required input:
* Input protein MSA directory
* Input CDS MSA directory
* Input protein exon position file
* Output protein MSA directory
* Output CDS MSA directory
* Output protein exon position file

### Exon position file format
The exon position file must have the following tab-delimited, one line per MSA:
* Name of MSA file (without directory) 
* The zero-based starting index position of each exon (protein coordinates)
* The length of the MSA (protein coordinates)

See the example given in the example folder for reference.

### MSA requirements
* Some non-trivial number of alignments must contain sequences from all species.
* The record id of each sequence (the defline until the first space) must be the species name.
* The file names for corresponding CDS and protein MSAs must be the same.
