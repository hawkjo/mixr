# MIXR: Mismatching Isoform eXon Remover

Required input:
* Input protein MSA directory
* Input CDS MSA directory
* Output protein MSA directory
* Output CDS MSA directory
* Protein exon position file

### Exon position file format
The exon position file must have the following tab-delimited, one line per MSA:
* Name of MSA file (without directory) 
* The zero-based starting index position of each exon
* The length of the MSA

See the example given in the example folder for reference.

### MSA requirements
* Some non-trivial number of alignments must contain sequences from all species.
* The record id of each sequence (the defline until the first space) must be the species name.
