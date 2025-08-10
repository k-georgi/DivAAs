# DivAAs (Divergent Amino Acids)
## Overview

DivAAs is a Python tool for processing, aligning, and comparing sequence groups provided in FASTA format. It supports nucleotide and peptide sequences, translating nucleotide input into amino acids if needed. The tool performs multiple sequence alignments using MAFFT, Clustal Omega, or MUSCLE, and computes conservation and mutation scores between groups. Scoring is based on substitution matrices (e.g., BLOSUM62) and hydrophobicity values. Results are saved to a CSV file and can be visualized with bar plots and sequence logo plots.

## Features

- **Sequence Processing**: Reads sequences from one or more FASTA files. Nucleotide sequences can be translated to amino acids using a specified genetic code range.
- **Group Assignment**: Assigns sequences to groups based on input files or an annotation file. Each input FASTA file can define a group, or you can provide 
                        a separate annotation (.txt) file to label sequences by group. Custom group names can be specified with the --group_names option.
- **Annotation Support**: Accepts annotation files from tools like ITOL or FigTree, or a custom text file with sequence IDs and group labels, to define group
                          membership.
- **Translation**: Translates nucleotide sequences to peptide sequences if requested, using provided start and stop positions.
- **Alignment**: Performs multiple sequence alignment using MAFFT, Clustal Omega, or MUSCLE (requires at least one tool installed and in your PATH).
- **Scoring**: Calculates conservation and mutation scores within and between groups, using BLOSUM62 scores and hydrophobicity for amino acids.
- **Visualization**: Generates bar plots (Matplotlib or Plotly) to display conservation and mutation scores. Optionally creates sequence logo plots for 
                     amino acids (--logo_plot y) and separate plots for nucleotide triplets (--nuc_plot y) or point mutations (--pm_plot y).
- **Output**: Saves results in a CSV file and can display summary plots. The CSV contains positional scores, most common residues or triplets per group,
              conservation rates, and mutation flags.

## Installation

### Prerequisites

- Python 3.x
- Required Python package: numpy
- Optional Python packages for graphics: matplotlib, plotly
- Optional for logo plots: tools like weblogo or logomaker (if using --logo_plot y)
- External alignment tools (at least one must be installed and in PATH): MAFFT, Clustal Omega, or MUSCLE.

### Installing Required Packages

Install the Python dependencies using pip:

```bash
pip install numpy matplotlib plotly
```

(For sequence logo plots, you may also need to install a logo generation library.)

### Installing External Alignment Tools

Ensure that at least one of the following alignment tools is installed:

- **MAFFT**: [Installation Guide](https://mafft.cbrc.jp/alignment/software/)
- **Clustal Omega**: [Installation Guide](http://www.clustal.org/omega/)
- **MUSCLE**: [Installation Guide](https://www.drive5.com/muscle/)

## Usage

### Command-Line Arguments
DivAAs is run from the command line. Specify one or more FASTA files and any desired options:

```bash
python DivAAs.py <fasta1> [<fasta2> ...] [--annotation ANNOTATION] [--startstop y] [--start START] [--stop STOP] \
                 [--conservation_score y] [--ref_ids REF_IDs ...] [--top TOP] [--lb LB] [--ub UB] [--positions POSITIONS ...] \
                 [--refseq REFSEQ] [--save y] [--directory DIRECTORY] [--show y] [--use_plotly y] [--nuc_plot n] [--pm_plot n] \
                 [--logo_plot y] [--group_names GROUP_NAMES ...]
```
#### Required Arguments

- `fasta`: Paths to one or more FASTA files containing sequence groups (required).

#### Optional Arguments

- `--annotation`: Path to a separate annotation file (.txt) that defines group membership by sequence ID.
- `--startstop`: Translate sequences from a defined start to stop position (y or n, default n).
- `--start`: Translation start position (integer, default 0).
- `--stop`: Translation stop position (integer, default None).
- `--conservation_score`: Calculate scores with conservation (y or n, default y).
- `--ref_ids`: Reference sequence IDs for each group (space-separated list, one ID per group).
- `--top`: Number of top-scoring positions to display in output (integer, default 20).
- `--lb`: Lower bound filter for positions (integer, default None).
- `--ub`: Upper bound filter for positions (integer, default None).
- `--positions`: Specific sequence positions to display (space-separated integers).
- `--refseq`: Reference sequence group (set to 1 or 2, default 1).
- `--save`: Save results to file (y or n, default y).
- `--directory`: Directory to save results (default is the user's home directory).
- `--show`: Display results on screen (y or n, default n).
- `--use_plotly`: Display plots using Plotly (y or n, default y).
- `--nuc_plot`: Show separate plot for nucleotide triplets (y or n, default n).
- `--pm_plot`: Show separate plot for point mutations (y or n, default n).
- `--logo_plot`: Show amino acid sequence logo plot (y or n, default y).
- `--group_names`: Names for each group (space-separated, in the same order as input files).

### Annotation Files
If you want to define group membership via an annotation file, use the --annotation option with a text file that assigns each sequence to a group. The annotation file should list sequence IDs and their group labels. You can generate such a file by:
- iTOL: Export annotations via the "Colors and styles annotation" option.
- FigTree: Save annotations via File → Export Trees → Tree File Format: NEXUS (check "Include Annotations").
- Custom File: Create a .txt file where each line has the sequence ID in the first column and the group label in the last column.
Save the annotation as a .txt file and pass its path to --annotation.

### Example Usage

**Basic Usage with One or More FASTA Files**:
```bash
python DivAAs.py group1.fasta group2.fasta
```
(Processes two FASTA files as two groups.)

**Using an Annotation File**:
```bash
python DivAAs.py sequences.fasta --annotation annotations.txt
```
(Defines group labels via annotations.txt.)

**Specifying Translation Start and Stop**:
```bash
python DivAAs.py seqs.fasta --startstop y --start 10 --stop 100
```

**Specifying Reference IDs and Group Names**:
```bash
python DivAAs.py group1.fasta group2.fasta --ref_ids SeqID1 SeqID2 --group_names Control Treatment
```

**Saving Results to a Specific Directory**:
```bash
python DivAAs.py seqs.fasta --directory /path/to/results
```

**Displaying Results with Plotly**:
```bash
python DivAAs.py seqs.fasta --use_plotly y
```

**Displaying Nucleotide and Point Mutation Plots**:
```bash
python DivAAs.py seqs.fasta --nuc_plot y --pm_plot y
```

**Displaying Amino Acid Sequence Logo Plot**:
```bash
python DivAAs.py seqs.fasta --logo_plot y
```

## Output

- **CSV File**: DivAAs saves a CSV file (DivAAs_output.csv) in the specified directory (default is the home directory). The CSV includes:
- `Position`: Position in the multiple sequence alignment.
- `Positions_in_<REF>`: Corresponding position(s) in the reference sequence(s) (IDs from --ref_ids).
- `Score`: Normalized conservation/mutation score for the position.
- `Most_Common_<GroupName>`: Most common amino acid (or nucleotide triplet) in each group at that position.
- `Conservation_<GroupName>`: Conservation rate (fraction of sequences sharing the most common residue) in each group.
- `Most_Common_Triplet_<GroupName>`: (If nucleotide input) most common codon (triplet) in each group.
- `Point_Mutations`: Indicates whether differences between groups could be due to a single point mutation (Y/N).

- **Plots**: In addition to the CSV, DivAAs can generate bar plots of conservation and mutation scores, and sequence logo plots if requested. Plot display is       
             controlled by `--show` and `--use_plotly`.

## Troubleshooting
- **No Alignment Tool Found**: Ensure MAFFT, Clustal Omega, or MUSCLE is installed and in your system PATH. The tool requires at least one to perform alignments.
- **Missing Python Packages**: Install any missing packages via pip install numpy matplotlib plotly.
- **File Not Found**: Verify that paths to FASTA and annotation files are correct.
- **Invalid Arguments**: Double-check argument syntax. Use python DivAAs.py --help to see all options.# COPS_private