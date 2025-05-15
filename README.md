# CoSeD (Comparison of Sequence Differences) Tool

## Overview

CoSeD is a Python-based tool designed to process, align, and compare sequence groups provided in FASTA files. It supports both nucleotide and peptide sequences, performs multiple sequence alignment using tools like Clustal Omega, MUSCLE, or MAFFT, and calculates conservation and mutation scores between groups based on substitution matrices (e.g., BLOSUM62) and hydrophobicity. The results are output as a CSV file and can be visualized using bar plots.

## Features

- **Sequence Processing**: Reads and processes sequences from FASTA files.
- **Group Assignment**: Assigns sequences to groups based on headers, annotations, or separate files.
- **Translation**: Translates nucleotide sequences to peptide sequences using the genetic code.
- **Alignment**: Performs multiple sequence alignment using MAFFT, Clustal Omega, or MUSCLE.
- **Scoring**: Calculates conservation and mutation scores between groups using BLOSUM62 and hydrophobicity values.
- **Visualization**: Generates bar plots using Matplotlib or Plotly to visualize conservation and mutation scores.
- **Output**: Saves results as a CSV file and provides options for displaying plots.

## Installation

### Prerequisites

- Python 3.x
- Required Python packages: `numpy`
- Optional Python packages for graphic representation:`matplotlib`, `plotly`
- External alignment tools: MAFFT, Clustal Omega, or MUSCLE (at least one must be installed and accessible in the system PATH)

### Installing Required Packages

You can install the Python packages using `pip`:

```bash
pip install numpy matplotlib plotly
```

### Installing External Alignment Tools

Ensure that at least one of the following alignment tools is installed:

- **MAFFT**: [Installation Guide](https://mafft.cbrc.jp/alignment/software/)
- **Clustal Omega**: [Installation Guide](http://www.clustal.org/omega/)
- **MUSCLE**: [Installation Guide](https://www.drive5.com/muscle/)

## Usage

### Command-Line Arguments

CoSeD can be run from the command line with the following arguments:

```bash
python CoSeD.py <file1> [--file2 FILE2] 
                 [--startstop STARTSTOP] [--start START] [--stop STOP] 
                 [--conservation_score CONSERVATION_SCORE] [--id1 ID1] [--id2 ID2] 
                 [--top TOP] [--lb LB] [--ub UB] [--positions POSITIONS]
                 [--refseq REFSEQ] [--save SAVE] [--directory DIRECTORY] [--show DISPLAY] 
                 [--use_plotly USE_PLOTLY] [--nuc_plot NUC_PLOT] [--pm_plot PM_PLOT] 
                 [--name1 NAME1] [--name2 NAME2]
```

#### Required Arguments

- `file1`: Path to the first FASTA file.

#### Optional Arguments

- `--file2`: Path to the second file (FASTA or TXT).
- `--startstop`: Translate sequence from defined start to defined stop (y/n). Default is 'n'.
- `--start`: Optional translation start position. Default is 0.
- `--stop`: Optional translation stop position. Default is None.
- `--conservation_score`: Calculate score with conservation or not (y/n). Default is 'y'.
- `--id1`: ID for reference sequence 1. Default is None.
- `--id2`: ID for reference sequence 2. Default is None.
- `--top`: Number of top scores to display. Default is 20.
- `--lb`: Lower bound for filtering positions. Default is None.
- `--ub`: Upper bound for filtering positions. Default is None.
- `--positions`: Defined positions to display, separeted by space. Default is None.
- `--refseq`: Reference sequence (1 or 2). Default is '1'.
- `--save`: Save results? (y/n). Default is 'y'.
- `--directory`: Directory to save results. Default is the user's home directory.
- `--show`: Display results? (y/n). Default is 'y'.
- `--use_plotly`: Display results with Plotly? (y/n). Default is 'y'.
- `--nuc_plot`: Show separate graph with triplets (Plotly only)? (y/n). Default is 'n'.
- `--pm_plot`: Show separate graph with point mutations? (y/n). Default is 'n'.
- `--name1`: Name of the first group. Default is 'Group 1'.
- `--name2`: Name of the second group. Default is 'Group 2'.

### Annotation Files
If you want to distinguish your groups with an annotation file, please use one of the following options:
- ITOL: export annotations via the option "Colors and styles annotation"
- figtree: save annotation file via File -> Export Trees -> Tree File Format: NEXUS (check box "Include Annotations")
- Create an individual annotation file with sequences in the first and the group annotation in the last column
The annotation file needs to be saved as '.txt' file

### Example Usage

1. **Basic Usage with One FASTA File**:
   ```bash
   python CoSeD.py sequences.fasta
   ```

2. **Using a Second FASTA File**:
   ```bash
   python CoSeD.py sequences1.fasta --file2 sequences2.fasta
   ```

3. **Using an Annotation File**:
   ```bash
   python CoSeD.py sequences.fasta --file2 annotations.txt
   ```

4. **Specifying Translation Start and Stop**:
   ```bash
   python CoSeD.py sequences.fasta --startstop y --start 10 --stop 100
   ```

5. **Saving Results to a Specific Directory**:
   ```bash
   python CoSeD.py sequences.fasta --directory /path/to/save/results
   ```

6. **Displaying Results with Plotly**:
   ```bash
   python CoSeD.py sequences.fasta --use_plotly y
   ```

7. **Displaying Nucleotide and Point Mutation Plots**:
   ```bash
   python CoSeD.py sequences.fasta --nuc_plot y --pm_plot y
   ```

## Output

- **CSV File**: The results are saved as `CoSeD_output.csv` in the specified directory. The file contains the following columns:
  - `Position`: Position in the alignment.
  - `Positions_in_<ID>`: Positions in the reference sequence.
  - `Score`: Normalized score for the position.
  - `Most_Common_<Group>`: Most common amino acid in each group.
  - `Conservation_<Group>`: Conservation rate in each group.
  - `Most_Common_Triplet_<Group>`: Most common nucleotide triplet in each group (if nucleotide sequences are provided).
  - `Point_Mutations`: Indicates if the difference could be due to a point mutation (Y/N).

- **Plots**: The tool generates bar plots showing conservation and mutation scores.

## Troubleshooting

- **No Alignment Tool Found**: Ensure that at least one of the alignment tools (MAFFT, Clustal Omega, or MUSCLE) is installed and accessible in your system PATH.
- **Missing Python Packages**: Install the required Python packages using `pip install numpy matplotlib plotly`.
- **File Not Found**: Ensure that the paths to the input files are correct.