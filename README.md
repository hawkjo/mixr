# MIXR: Mismatching Isoform eXon Remover

A package for the removal of exons with mismatching isoforms as described in:

### A meta-analysis of bat phylogenetics based on whole genomes and transcriptomes of 18 species

**John A. Hawkins, Maria E. Kaczmarek, William H. Press, Sara L. Sawyer**

(in review)

### Installation

The following instructions should work across platforms, except that installing virtualenv with apt-get is Ubuntu specific. For other platforms, install virtualenv appropriately if desired.

First, clone the repository to a local directory:

```
git clone https://github.com/hawkjo/mixr.git
```

Optionally, you can install into a virtual environment (recommended):

```
sudo apt-get install -y virtualenv
cd mixr
virtualenv envmixr
. envmixr/bin/activate
```

Now install required packages and mixr using pip:

```
pip install -r requirements.txt && python setup.py install
```

### Usage

```
Usage:
  mixr <in_prot_dir> <in_cds_dir> <in_exon_pos_file> <out_prot_dir> <out_cds_dir> <out_exon_pos_file> [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.
```

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

See the example given in the examples folder.

### MSA requirements
* Some non-trivial number of alignments must contain sequences from all species.
* The record id of each sequence (the defline until the first space) must be the species name.
* The file names for corresponding CDS and protein MSAs must be the same.
