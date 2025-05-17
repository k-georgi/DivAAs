"""
Description:
--------------------
This script processes, aligns, and compares sequence groups provided in FASTA files. It reads
FASTA and annotation text files (from sources such as ITOL or FigTree), groups sequences,
performs multiple sequence alignment using available tools (Clustal Omega, MUSCLE, or MAFFT),
calculates conservation and mutation scores between groups based on substitution matrices
(e.g., BLOSUM62) and hydrophobicity, and finally outputs the results as a CSV file and a bar plot.
"""

# region Import
import argparse                     # For command-line argument parsing
import collections
import importlib.util               # For checking availability of Python packages
import logomaker as lm
import matplotlib.pyplot as plt     # For plotting results with matplotlib
import numpy as np                  # For numerical operations and array handling
import os                           # For file system operations
import pandas as pd
import plotly.graph_objects as go   # For plotting results with plotly
import random                       # For random selection in case of ties
import re                           # For regular expression matching
import shutil                       # For checking availability of external tools
import subprocess                   # For running external alignment tools
import tempfile                     # For creating temporary files
#from collections import Counter     # For counting elements in sequences
from typing import Dict, List, Tuple, Optional, Union, Any
# endregion

# region Constants
class Constants:
    GENETIC_CODE = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATN': 'I',  # Isoleucin
        'ATG': 'M',  # Methionin (Start-Codon)
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'ACN': 'T',  # Threonin
        'AAC': 'N', 'AAT': 'N',  # Asparagin
        'AAA': 'K', 'AAG': 'K',  # Lysin
        'AGC': 'S', 'AGT': 'S',  # Serin
        'AGA': 'R', 'AGG': 'R',  # Arginin
        
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CTN': 'L',  # Leucin
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CCN': 'P',  # Prolin
        'CAC': 'H', 'CAT': 'H',  # Histidin
        'CAA': 'Q', 'CAG': 'Q',  # Glutamin
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CGN': 'R',  # Arginin
        
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GTN': 'V',  # Valin
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GCN': 'A',  # Alanin
        'GAC': 'D', 'GAT': 'D',  # Asparaginsäure
        'GAA': 'E', 'GAG': 'E',  # Glutaminsäure
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GGN': 'G',  # Glycin
        
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TCN': 'S',  # Serin
        'TTC': 'F', 'TTT': 'F',  # Phenylalanin
        'TTA': 'L', 'TTG': 'L',  # Leucin
        'TAC': 'Y', 'TAT': 'Y',  # Tyrosin
        'TAA': '*', 'TAG': '*', 'TGA': '*',  # Stop-Codons
        'TGC': 'C', 'TGT': 'C',  # Cystein
        'TGG': 'W'  # Tryptophan
    }

    STOP_CODONS = ["TAG", "TGA", "TAA"]

    START_CODON = 'ATG'

    BLOSUM62 = {
        ('A', 'A'): 4, ('A', 'C'): 0, ('A', 'D'): -2, ('A', 'E'): -1, ('A', 'F'): -2,
        ('A', 'G'): 0, ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'K'): -1, ('A', 'L'): -1,
        ('A', 'M'): -1, ('A', 'N'): -2, ('A', 'P'): -1, ('A', 'Q'): -1, ('A', 'R'): -1,
        ('A', 'S'): 1, ('A', 'T'): 0, ('A', 'V'): 0, ('A', 'W'): -3, ('A', 'Y'): -2,
        
        ('C', 'C'): 9, ('C', 'D'): -3, ('C', 'E'): -4, ('C', 'F'): -2, ('C', 'G'): -3,
        ('C', 'H'): -3, ('C', 'I'): -1, ('C', 'K'): -3, ('C', 'L'): -1, ('C', 'M'): -1,
        ('C', 'N'): -3, ('C', 'P'): -3, ('C', 'Q'): -3, ('C', 'R'): -3, ('C', 'S'): -1,
        ('C', 'T'): -1, ('C', 'V'): -1, ('C', 'W'): -2, ('C', 'Y'): -2,

        ('D', 'D'): 6, ('D', 'E'): 2, ('D', 'F'): -3, ('D', 'G'): -1, ('D', 'H'): -1,
        ('D', 'I'): -3, ('D', 'K'): -1, ('D', 'L'): -4, ('D', 'M'): -3, ('D', 'N'): 1,
        ('D', 'P'): -1, ('D', 'Q'): 0, ('D', 'R'): -2, ('D', 'S'): 0, ('D', 'T'): -1,
        ('D', 'V'): -3, ('D', 'W'): -4, ('D', 'Y'): -3,

        ('E', 'E'): 5, ('E', 'F'): -3, ('E', 'G'): -2, ('E', 'H'): 0, ('E', 'I'): -3,
        ('E', 'K'): 1, ('E', 'L'): -3, ('E', 'M'): -2, ('E', 'N'): 0, ('E', 'P'): -1,
        ('E', 'Q'): 2, ('E', 'R'): 0, ('E', 'S'): 0, ('E', 'T'): -1, ('E', 'V'): -2,
        ('E', 'W'): -3, ('E', 'Y'): -2,

        ('F', 'F'): 6, ('F', 'G'): -3, ('F', 'H'): -1, ('F', 'I'): 0, ('F', 'K'): -3,
        ('F', 'L'): 0, ('F', 'M'): 0, ('F', 'N'): -3, ('F', 'P'): -4, ('F', 'Q'): -3,
        ('F', 'R'): -3, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'V'): -1, ('F', 'W'): 1,
        ('F', 'Y'): 3,
    }

    # Hydrophobicity values for amino acids (Kyle and Doolittle scale)
    hydrophobicity = {
        'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
        'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
        'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
        'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
    }
# endregion

# region FASTA and Sequence-Type
class FASTAProcessor:
    # Define constants at the beginning of the class
    NUCLEOTIDES: frozenset[str] = frozenset('ATGCNUatgcnu')
    AMINO_ACIDS: frozenset[str] = frozenset('DEFHIKLMPQRSVWYdefhiklmpqrsvwy')

    
    # Sequence type constants
    TYPE_NUCLEOTIDE: str = 'NUCLEOTIDE'
    TYPE_PEPTIDE: str = 'PEPTIDE'
    TYPE_UNKNOWN: str = 'UNKNOWN'

    
    def __init__(self, filepath: str):
        """
        Initializes a FASTA file object.
        
        Args:
            filepath (str): Path to the FASTA file
        """
        self.filepath = filepath
        self.sequences: Dict[str, str] = {}
        self._load_sequences()
    
    def _load_sequences(self):
        """
        Loads sequences from the FASTA file.
        """
        try:
            with open(self.filepath, 'r') as file:
                current_header: str = ""
                current_sequence: List[str] = [] # Use a list for efficient string concatenation
                
                for line in file:
                    line = line.strip()
                    if line.startswith('>'):
                        # Store the previous sequence if it exists
                        if current_header and current_sequence:
                            self.sequences[current_header] = ''.join(current_sequence)
                        current_header = line[1:]  # Remove '>' from the header
                        current_sequence = []
                    elif line:  # Only add non-empty lines to the sequence
                        current_sequence.append(line)
                
                # Store the last sequence
                if current_header and current_sequence:
                    self.sequences[current_header] = ''.join(current_sequence)
                    
        except FileNotFoundError:
            raise FileNotFoundError(f"The file {self.filepath} was not found.")
    
    def _determine_type(self, sequence: str) -> str:
        """
        Determines whether the sequence is a nucleotide or peptide sequence.
        
        Args:
            sequence (str): The sequence to analyze
            
        Returns:
            str: The type of the sequence (NUCLEOTIDE, PEPTIDE, or UNKNOWN)
        """
        # Clean sequence: remove spaces, tabs, and newlines
        clean_seq = ''.join(c for c in sequence if c not in {' ', '\t', '\n'})
        
        if not clean_seq:
            return self.TYPE_UNKNOWN
        
        # Count the number of characters
        total_chars = len(clean_seq)
        nucleotide_chars = sum(1 for c in clean_seq if c in self.NUCLEOTIDES)
        
        # Decision criterion: If more than 90% of the characters are nucleotides, classify as nucleotide sequence
        if nucleotide_chars / total_chars > 0.9:
            return self.TYPE_NUCLEOTIDE
        else:
            return self.TYPE_PEPTIDE
    
    @property
    def type(self):
        """
        Determines the dominant sequence type in the file.
        
        Returns:
            str: The dominant sequence type (NUCLEOTIDE, PEPTIDE, or UNKNOWN)
        """
        if not self.sequences:
            return self.TYPE_UNKNOWN
        
        # Count occurrences of each sequence type
        type_counts: Dict[str, int] = {self.TYPE_NUCLEOTIDE: 0, self.TYPE_PEPTIDE: 0, self.TYPE_UNKNOWN: 0}
        
        for seq in self.sequences.values():
            seq_type = self._determine_type(seq)
            type_counts[seq_type] += 1
        
        # Determine the dominant type
        if type_counts[self.TYPE_NUCLEOTIDE] > type_counts[self.TYPE_PEPTIDE]:
            return self.TYPE_NUCLEOTIDE
        elif type_counts[self.TYPE_PEPTIDE] > type_counts[self.TYPE_NUCLEOTIDE]:
            return self.TYPE_PEPTIDE
        else:
            return self.TYPE_UNKNOWN

# endregion

# region Assign Groups
class GroupAssigner:
    """
    Assigns sequences to group A or group B based on different criteria.
    """
    def __init__(self):
        """Initialize empty dictionaries for group A and group B."""
        self.group_A: Dict[str, str] = {}
        self.group_B: Dict[str, str] = {}
    
    def assign_by_headers(self, sequences: Dict[str, str]) -> Optional[Tuple[Dict[str, str], Dict[str, str]]]:
        """
        Assigns sequences to groups based on the last character of their headers.
        
        Args:
            sequences: Dictionary mapping sequence headers to sequence data
            
        Returns:
            Tuple of two dictionaries (group_A, group_B) with sequences assigned to each group,
            or None if sequences dictionary is empty
        """
        if not sequences:
            return None
        
        # Reset groups to ensure clean assignment
        self.group_A.clear()
        self.group_B.clear()
        
        # Use the last character of the first header as reference for group A
        headers = list(sequences.keys())
        first_group = headers[0].strip()[-1]
        
        for header, seq in sequences.items():
            clean_header = header.strip()
            base_header = clean_header.split()[0]
            
            if clean_header[-1] == first_group:
                self.group_A[f"{base_header} A"] = seq
            else:
                self.group_B[f"{base_header} B"] = seq
        
        return self.group_A, self.group_B
    
    def assign_by_annotation(self, 
                            sequences: Dict[str, str], 
                            annotations: Dict[str, str]) -> Optional[Tuple[Dict[str, str], Dict[str, str]]]:
        """
        Assigns sequences to groups based on annotation data.
        
        Args:
            sequences: Dictionary mapping sequence headers to sequence data
            annotations: Dictionary mapping sequence IDs to group annotations
            
        Returns:
            Tuple of two dictionaries (group_A, group_B) with sequences assigned to each group,
            or None if either input dictionary is empty
        """
        if not sequences or not annotations:
            return None
        
        # Reset groups to ensure clean assignment
        self.group_A.clear()
        self.group_B.clear()
        
        # Get first annotation as reference for group A
        first_group = None
        
        for header, seq in sequences.items():
            # Extract sequence ID from header (first word)
            seq_id = header.split()[0]
            
            # Set reference group if not already set
            if first_group is None and seq_id in annotations:
                first_group = annotations[seq_id]
            
            if seq_id in annotations:
                if annotations[seq_id] == first_group:
                    self.group_A[f"{seq_id} A"] = seq
                else:
                    self.group_B[f"{seq_id} B"] = seq
            else:
                # Default to group A if no annotation found
                self.group_A[f"{seq_id} A"] = seq
        
        return self.group_A, self.group_B
    
    def assign_separate_files(self, 
                             file1_sequences: Dict[str, str], 
                             file2_sequences: Dict[str, str]) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Assigns sequences from two separate files to groups A and B respectively.
        
        Args:
            file1_sequences: Dictionary of sequences from first file, assigned to group A
            file2_sequences: Dictionary of sequences from second file, assigned to group B
            
        Returns:
            Tuple of two dictionaries (group_A, group_B) with sequences assigned to each group
        """
        # Reset groups to ensure clean assignment
        self.group_A.clear()
        self.group_B.clear()
        
        # Add all sequences from file1 to group A
        for header, seq in file1_sequences.items():
            base_header = header.split()[0]
            self.group_A[f"{base_header} A"] = seq
        
        # Add all sequences from file2 to group B
        for header, seq in file2_sequences.items():
            base_header = header.split()[0]
            self.group_B[f"{base_header} B"] = seq
        
        return self.group_A, self.group_B


class AnnotationReader:
    """
    Reads and parses annotation files in various formats.
    """
    def __init__(self, txt_file: str):
        """
        Initializes the AnnotationReader with a text file path.
        
        Args:
            txt_file: Path to the annotation file
        """
        self.filepath = txt_file
        self.annotation = self.read_txt()
    
    def read_txt(self) -> Dict[str, str]:
        """
        Loads annotations from the file and automatically determines the file type.
        
        Returns:
            Dictionary with sequence IDs as keys and group annotations as values
        """
        try:
            with open(self.filepath, "r") as file:
                content = file.read()
                
                # Determine file type based on keywords
                if "DATA" in content.upper():
                    return self._read_itol()
                elif "NEXUS" in content.upper():
                    return self._read_figtree()
                elif "ANNOTATION" in content.upper():
                    return self._read_annotation()
                
                # If no keyword found, try as simple annotation file
                return self._read_annotation()
        except Exception as e:
            print(f"Error reading annotation file {self.filepath}: {e}")
            return {}

    def _read_itol(self) -> Dict[str, str]:
        """
        Reads an iTOL annotation file and returns a dictionary.
        
        Returns:
            Dictionary with sequence IDs as keys and groups as values
        """
        annotations: Dict[str, str] = {}
        data_started = False
        
        try:
            with open(self.filepath, "r") as file:
                for line in file:
                    line = line.strip()
                    if not line:
                        continue  # Skip empty lines
                        
                    if not data_started:
                        if "DATA" in line.upper():
                            data_started = True
                        continue  # Skip header lines until data begins
                        
                    # Assume line contains space-separated fields
                    parts = line.split()
                    if len(parts) >= 2:
                        seq_id = parts[0]     # First element is the sequence ID
                        group = parts[-1]     # Last element should be the group
                        annotations[seq_id] = group
        except Exception as e:
            print(f"Error reading iTOL file {self.filepath}: {e}")
            
        return annotations

    def _read_figtree(self) -> Dict[str, str]:
        """
        Reads a FigTree annotation file and extracts identifiers with group assignments.
        
        Returns:
            Dictionary with sequence IDs as keys and groups as values
        """
        try:
            with open(self.filepath, "r") as f:
                data = f.read()

            # Try to extract annotations with !name
            pattern_name = r"(?<!\))'([^']+)'\[.*?!name=\"(.*?)\""
            matches = re.findall(pattern_name, data)

            if matches:
                # Create dictionary from matches: identifier as key, group as value
                return {identifier: group for identifier, group in matches}

            # If no !name annotations found, try with !color
            pattern_color = r"(?<!\))'([^']+)'\[.*?!color=#([0-9a-fA-F]{6})"
            matches_color = re.findall(pattern_color, data)

            if matches_color:
                annotations: Dict[str, str] = {}
                unique_colors: List[str] = []
                
                for identifier, color in matches_color:
                    # Register unique color if not seen before
                    if color not in unique_colors:
                        unique_colors.append(color)
                    idx = unique_colors.index(color)
                    # Assign group letter based on order (0 -> A, 1 -> B, etc.)
                    group = chr(65 + idx)  # 65 corresponds to "A" in ASCII
                    annotations[identifier] = group
                    
                return annotations

            # If no annotations found, return empty dictionary
            return {}
        except Exception as e:
            print(f"Error reading FigTree file {self.filepath}: {e}")
            return {}

    def _read_annotation(self) -> Dict[str, str]:
        """
        Reads a simple annotation file with identifiers and group assignments.
        
        Returns:
            Dictionary with identifiers and group assignments
        """
        annotations: Dict[str, str] = {}
        try:
            with open(self.filepath, "r") as file:
                for line in file:
                    line = line.strip()
                    if not line or line.upper().startswith("ANNOTATION"):
                        continue  # Skip empty lines and header
                        
                    parts = line.split()
                    if len(parts) >= 2:
                        seq_id = parts[0]     # First part is assumed to be the sequence ID
                        group = parts[-1]     # Last part should be the group
                        annotations[seq_id] = group
        except Exception as e:
            print(f"Error reading annotation file {self.filepath}: {e}")
            
        return annotations

class FileProcessor:
    """
    Processes FASTA files and assigns sequences to groups based on different criteria.
    """
    def __init__(self):
        """Initialize with a GroupAssigner."""
        self.assigner = GroupAssigner()
        self.fasta_file1: Dict[str, str] = {}
    
    def process_files(self, 
                     file1: str, 
                     file2_or_filetxt: Optional[str] = None) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Processes at least one FASTA file (file1) and optionally a second file.
        
        Args:
            file1: Path to the first FASTA file (required)
            file2_or_filetxt: Path to a second file, either FASTA or TXT (optional)
        
        Returns:
            Two dictionaries (group_A, group_B) with the assigned sequences
            
        Raises:
            ValueError: If file1 is not provided or if file2_or_filetxt has an invalid format
        """
        if not file1:
            raise ValueError("Error: A FASTA file (file1) must be specified!")
        
        # Load the first FASTA file
        self.fasta_file1 = FASTAProcessor(file1)
        
        # Determine type and behavior based on the second file
        if file2_or_filetxt:
            if file2_or_filetxt.lower().endswith((".fasta", ".fa", ".fas")):
                # Fall 1: Zweite FASTA-Datei
                fasta_file2 = FASTAProcessor(file2_or_filetxt)
                return self.assigner.assign_separate_files(self.fasta_file1.sequences, fasta_file2.sequences)
            
            elif file2_or_filetxt.lower().endswith(".txt"):
                # Fall 2a: Annotationsdatei
                annotations = AnnotationReader(file2_or_filetxt)
                return self.assigner.assign_by_annotation(self.fasta_file1.sequences, annotations.annotation)
            
            else:
                raise ValueError("Error: Die zweite Datei muss entweder eine FASTA-Datei (.fasta, .fa, .fas) oder eine TXT-Datei (.txt) sein.")
        
        else:
            # Case 2b: Group assignment based on headers
            result = self.assigner.assign_by_headers(self.fasta_file1.sequences)
            if result is None:
                return {}, {}
            return result

# endregion

# region Translation
class TranslateSequence:
    """
    Translates nucleotide sequences to peptide sequences using different methods.
    
    This class provides functionality to translate DNA sequences to protein sequences
    either by finding start codons and reading frames automatically or by specifying
    exact start and end positions for translation.
    """
    
    def __init__(self, 
                seq_dict: Dict[str, str], 
                StartStop: str = 'N', 
                start: Optional[int] = 0, 
                stop: Optional[int] = None):
        """
        Initialize the sequence translator.
        
        Args:
            seq_dict: Dictionary of sequence headers and their nucleotide sequences
            StartStop: 'Y'/'YES' to use explicit start/stop positions, 'N' to find open reading frames
            start: Starting position for translation when StartStop is enabled
            stop: Ending position for translation when StartStop is enabled (None for end of sequence)
        """
        self.GENETIC_CODE = Constants.GENETIC_CODE
        self.STOP_CODONS = Constants.STOP_CODONS
        self.START_CODON = Constants.START_CODON
        
        self.seq_dict = seq_dict
        self.StartStop = StartStop.upper()
        self.start = start
        self.stop = stop
        
        self.nuc_seq = {}
        self.pep_seq = {}
        
        self._process_sequences()
    
    def _process_sequences(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Process all sequences in seq_dict and generate nucleotide and peptide dictionaries.
        
        Returns:
            A tuple containing (nucleotide_sequences, peptide_sequences) dictionaries
        """
        use_startstop = self.StartStop in ['Y', 'YES']
        
        for header, seq in self.seq_dict.items():
            try:
                # Select translation method based on StartStop parameter
                if use_startstop:
                    nuc, pep = self._translate_startstop(seq, self.start, self.stop)
                else:
                    nuc, pep = self._translate(seq)
                
                # Only store sequences that were successfully translated
                if nuc and pep:
                    self.nuc_seq[header.split()[0]] = nuc
                    self.pep_seq[header] = pep
                    
            except Exception as e:
                print(f"Error translating sequence '{header}': {e}")

        return self.nuc_seq, self.pep_seq
    
    def _translate(self, seq: str) -> Tuple[str, str]:
        """
        Translate a nucleotide sequence to a peptide by finding the longest open reading frame.
        
        Searches for ATG start codons and translates until finding a stop codon,
        returning the longest complete translated sequence.
        
        Args:
            seq: Nucleotide sequence to translate
            
        Returns:
            Tuple of (nucleotide_sequence, peptide_sequence) for the longest ORF found
        """
        if not seq:
            return "", ""
        
        seq = seq.upper()
        
        longest_nuc = ""
        longest_pep = ""
        
        # Find all start codon positions
        start_positions = self._find_all_occurrences(seq, self.START_CODON)
        
        # Try translation from each start position
        for start_pos in start_positions:
            nuc_seq = ""
            pep_seq = ""
            pos = start_pos
            
            # Translate codon by codon until reaching a stop codon or end of sequence
            while pos <= len(seq) - 3:
                codon = seq[pos:pos+3]
                
                # Handle unknown codons
                if codon in self.GENETIC_CODE:
                    aa = self.GENETIC_CODE[codon]
                else:
                    aa = "X"  # Standard symbol for unknown amino acid
                
                # Stop at stop codon
                if aa == "*":
                    break
                
                nuc_seq += codon
                pep_seq += aa
                pos += 3
            
            # Keep track of the longest translation
            if len(nuc_seq) > len(longest_nuc):
                longest_nuc = nuc_seq
                longest_pep = pep_seq
        
        return longest_nuc, longest_pep
    
    def _translate_startstop(self, seq: str, start: int, stop: Optional[int] = None) -> Tuple[str, str]:
        """
        Translate a nucleotide sequence using specific start and stop positions.
        
        Args:
            seq: Nucleotide sequence to translate
            start: Starting position for translation (0-based index)
            stop: Ending position for translation (None for end of sequence)
            
        Returns:
            Tuple of (nucleotide_sequence, peptide_sequence)
        """
        if not seq:
            return "", ""
        
        seq = seq.upper()
        nuc_seq = ""
        pep_seq = ""
        
        # Set stop to end of sequence if not specified
        if stop is None:
            stop = len(seq)
        
        # Ensure we only process up to the specified stop position
        effective_stop = min(stop, len(seq))
        
        # Translate codon by codon
        pos = start
        while pos <= effective_stop - 3:
            codon = seq[pos:pos+3]
            
            # Handle unknown codons
            if codon in self.GENETIC_CODE:
                aa = self.GENETIC_CODE[codon]
            else:
                aa = "X"  # Standard symbol for unknown amino acid
            
            # Stop at stop codon
            if aa == "*":
                break
            
            nuc_seq += codon
            pep_seq += aa
            pos += 3
        
        return nuc_seq, pep_seq
    
    def _find_all_occurrences(self, seq: str, pattern: str) -> List[int]:
        """
        Find all occurrences of a pattern in a sequence.
        
        Args:
            seq: The sequence to search in
            pattern: The pattern to search for
            
        Returns:
            List of starting positions where the pattern was found
        """
        positions = []
        start = 0
        
        while True:
            start = seq.find(pattern, start)
            if start == -1:
                break
            positions.append(start)
            start += 1
            
        return positions

# endregion
        
# region Alignment
class Alignment:
    """
    A class for performing and managing pairwise sequence alignments.
    
    This class handles the alignment of protein sequences with optional corresponding
    nucleotide sequences. It supports multiple alignment tools (MAFFT, Clustal Omega, MUSCLE)
    and provides functionality to map nucleotide sequences to protein alignments.
    
    Attributes:
        pep_A (Dict[str, str]): Dictionary of protein sequences for group A
        pep_B (Dict[str, str]): Dictionary of protein sequences for group B
        nuc_A (Dict[str, str], optional): Dictionary of nucleotide sequences for group A
        nuc_B (Dict[str, str], optional): Dictionary of nucleotide sequences for group B
        pep_aligned_A (Dict[str, str]): Dictionary of aligned protein sequences for group A
        pep_aligned_B (Dict[str, str]): Dictionary of aligned protein sequences for group B
        nuc_aligned_A (Dict[str, str], optional): Dictionary of nucleotide sequences mapped to protein alignments for group A
        nuc_aligned_B (Dict[str, str], optional): Dictionary of nucleotide sequences mapped to protein alignments for group B
    """
    
    def __init__(self, pep_A: Dict[str, str], pep_B: Dict[str, str], 
                 nuc_A: Optional[Dict[str, str]] = None, nuc_B: Optional[Dict[str, str]] = None):
        """
        Initialize an Alignment object and perform the alignment.
        
        Args:
            pep_A: Dictionary of protein sequences for group A
            pep_B: Dictionary of protein sequences for group B
            nuc_A: Dictionary of nucleotide sequences for group A (optional)
            nuc_B: Dictionary of nucleotide sequences for group B (optional)
        """
        self.pep_A = pep_A
        self.pep_B = pep_B
        self.nuc_A = nuc_A
        self.nuc_B = nuc_B
        self.pep_aligned_A: Dict[str, str] = {}
        self.pep_aligned_B: Dict[str, str] = {}
        self.nuc_aligned_A: Optional[Dict[str, str]] = {}
        self.nuc_aligned_B: Optional[Dict[str, str]] = {}
        
        # Perform protein alignment
        self._run_alignment()
        
        # Map nucleotide sequences to protein alignment if provided
        if nuc_A and nuc_B:
            self.nuc_aligned_A = self._map_nucleotides_to_alignment(self.nuc_A, self.pep_aligned_A)
            self.nuc_aligned_B = self._map_nucleotides_to_alignment(self.nuc_B, self.pep_aligned_B)

    def _merge_seq(self, sequences_A: Dict[str, str], sequences_B: Dict[str, str]) -> Dict[str, Dict[str, str]]:
        """
        Merge two dictionaries containing FASTA sequences from groups A and B.

        Args:
            sequences_A: Dictionary with the first group of sequences
            sequences_B: Dictionary with the second group of sequences

        Returns:
            A dictionary with keys "A" and "B" and corresponding sequence dictionaries as values.
        """
        return {
            "A": sequences_A,  # All sequences from group A
            "B": sequences_B   # All sequences from group B
        }
    
    @staticmethod
    def check_alignment_tools() -> Dict[str, bool]:
        """
        Check whether Clustal Omega, MUSCLE, or MAFFT is installed on the system.
        
        Returns:
            A dictionary with tool names as keys and booleans as values indicating if each tool is found.
        """
        tools = ["clustalo", "muscle", "mafft"]
        return {tool: shutil.which(tool) is not None for tool in tools}

    def _run_alignment(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Perform a multiple sequence alignment using an available alignment tool.
        
        This method checks for available alignment tools (MAFFT, Clustal Omega, or MUSCLE)
        and uses the first available one to perform the alignment.
        
        Returns:
            A tuple containing two dictionaries: (aligned_sequences_A, aligned_sequences_B)
        
        Raises:
            RuntimeError: If no alignment tool is installed.
        """
        # Merge Input Sequence Groups
        merged_dict = self._merge_seq(self.pep_A, self.pep_B)

        # Check Alignment Tools
        tools = self.check_alignment_tools()
        
        if not any(tools.values()):
            raise RuntimeError("No alignment tool installed. Please install MAFFT, Clustal Omega, or MUSCLE.")

        # Create input sequence records with a group suffix in the description
        input_records = []
        for group_key in ['A', 'B']:
            group_seqs = merged_dict.get(group_key, {})
            for seq_id, seq in group_seqs.items():
                new_desc = f"Seq_{seq_id}_{group_key}"  # Append group identifier to description
                input_records.append(f">{seq_id} {new_desc}\n{seq}\n")
        
        # Create temporary input and output files for the alignment tool
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as input_file:
            input_filename = input_file.name
            input_file.write("".join(input_records))
        
        output_filename = tempfile.mktemp(suffix='.fasta')
        
        try:
            self._execute_alignment_tool(tools, input_filename, output_filename)
            
            # Process the output file to parse aligned sequences
            aligned_sequences = self._parse_alignment_output(output_filename)
            
        except Exception as e:
            raise RuntimeError(f"Alignment failed: {e}")
        finally:
            # Clean up temporary files
            if os.path.exists(input_filename):
                os.remove(input_filename)
            if os.path.exists(output_filename):
                os.remove(output_filename)

        # Split the aligned sequences into group A and group B
        self.pep_aligned_A, self.pep_aligned_B = self._split_seq(aligned_sequences)
        
        return self.pep_aligned_A, self.pep_aligned_B
    
    def _execute_alignment_tool(self, tools: Dict[str, bool], input_filename: str, output_filename: str) -> None:
        """
        Execute the appropriate alignment tool based on what's available.
        
        Args:
            tools: Dictionary of available alignment tools
            input_filename: Path to the input file
            output_filename: Path to the output file
        
        Raises:
            RuntimeError: If execution of the alignment tool fails
        """
        if tools["mafft"]:
            # Set up and execute MAFFT command
            cmd = [
                "mafft",         # MAFFT command
                "--auto",        # Automatic alignment settings
                "--quiet",       # Suppress output
                input_filename   # Input file for alignment
            ]
            with open(output_filename, "w") as outfile:
                subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)
        elif tools["clustalo"]:
            # Set up and execute Clustal Omega command
            cmd = [
                "clustalo",           # Clustal Omega command
                "-i", input_filename,  # Input file
                "-o", output_filename, # Output file
                "--outfmt=fasta",      # Output format: FASTA
                "--force",             # Force execution even if output exists
                "--auto"               # Automatic settings
            ]
            subprocess.run(cmd, check=True)
        elif tools["muscle"]:
            # Set up and execute MUSCLE command
            cmd = [
                "muscle",
                "-align", input_filename,
                "-output", output_filename,
                "-quiet"
            ]
            subprocess.run(cmd, check=True)
    
    def _parse_alignment_output(self, output_filename: str) -> Dict[str, Tuple[str, str]]:
        """
        Parse the alignment output file and extract aligned sequences.
        
        Args:
            output_filename: Path to the alignment output file
        
        Returns:
            Dictionary mapping sequence IDs to tuples of (description, aligned sequence)
        """
        aligned_sequences = {}
        with open(output_filename, 'r') as output_file:
            current_id = None
            current_desc = None
            current_seq = []
            for line in output_file:
                line = line.strip()
                if line.startswith(">"):  # Header line
                    if current_id is not None:
                        aligned_sequences[current_id] = (current_desc, "".join(current_seq))
                    header_parts = line[1:].split(maxsplit=1)
                    current_id = header_parts[0]
                    current_desc = header_parts[1] if len(header_parts) > 1 else ""
                    current_seq = []
                else:
                    current_seq.append(line)
            # Don't forget the last sequence
            if current_id is not None:
                aligned_sequences[current_id] = (current_desc, "".join(current_seq))
        
        return aligned_sequences
    
    def _split_seq(self, sequences: Dict[str, Tuple[str, str]]) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Split the input aligned sequences into two dictionaries based on group identifiers in descriptions.

        Args:
            sequences: Dictionary with aligned sequences

        Returns:
            Two dictionaries: one for group "A" and one for group "B".
        """
        group_A: Dict[str, str] = {}  # Dictionary for group A sequences
        group_B: Dict[str, str] = {}  # Dictionary for group B sequences

        for seq_id, (description, sequence) in sequences.items():
            group_key = description.strip()[-1]  # Get the last character which indicates the group
            seq_id_clean = seq_id.split()[0]  # Remove any description part from the ID
            
            if group_key == "A":
                group_A[seq_id_clean] = sequence
            elif group_key == "B":
                group_B[seq_id_clean] = sequence
                
        return group_A, group_B

    def _map_nucleotides_to_alignment(self, nuc_sequences: Dict[str, str], pep_alignment: Dict[str, str]) -> Dict[str, str]:
        """
        Map original nucleotide triplets to aligned amino acids.
        
        This method creates a nucleotide alignment based on the protein alignment by
        expanding each amino acid position to its corresponding codon (3 nucleotides).
        Gaps in the protein alignment are expanded to 3 nucleotide gaps.
        
        Args:
            nuc_sequences: Original nucleotide sequences
            pep_alignment: Aligned protein sequences
        
        Returns:
            Dictionary with sequence IDs and corresponding nucleotide alignments
        """
        nuc_alignment: Dict[str, str] = {}

        for header in pep_alignment:
            pep_seq = pep_alignment[header]
            nuc_seq = nuc_sequences[header]
            
            # Extract nucleotide triplets
            nuc_mapping: List[str] = []
            nuc_pos: int = 0
            
            for aa in pep_seq:
                if aa == '-':
                    # For gaps in protein sequence, add three gap characters
                    nuc_mapping.append('-' * 3)
                else:
                    # For amino acids, get the corresponding nucleotide triplet
                    triplet = nuc_seq[nuc_pos:nuc_pos+3]
                    nuc_mapping.append(triplet)
                    nuc_pos += 3
            
            nuc_alignment[header] = ''.join(nuc_mapping)
        
        return nuc_alignment
          
# endregion
    
# region Score
class Scoring:
    """
    A class for calculating similarity scores between two groups of aligned sequences.
    
    This class compares aligned protein sequences from two groups (A and B), with optional
    nucleotide sequences, and calculates position-specific scores based on amino acid
    differences, BLOSUM62 substitution matrix, hydrophobicity, and conservation rates.
    
    Attributes:
        BLOSUM62 (Dict): BLOSUM62 substitution matrix for scoring amino acid exchanges
        hydrophobicity (Dict): Dictionary mapping amino acids to their hydrophobicity values
        pep_A (Dict[str, str]): Dictionary of aligned protein sequences for group A
        pep_B (Dict[str, str]): Dictionary of aligned protein sequences for group B
        nuc_A (Optional[Dict[str, str]]): Dictionary of aligned nucleotide sequences for group A
        nuc_B (Optional[Dict[str, str]]): Dictionary of aligned nucleotide sequences for group B
        idA (Optional[str]): ID of the reference sequence for group A
        idB (Optional[str]): ID of the reference sequence for group B
        con_score (str): Flag indicating whether to include conservation in score calculation ("y"/"n")
        results (np.ndarray): Array containing the calculated comparison results
    """
    
    def __init__(
        self, 
        pep_A: Dict[str, str], 
        pep_B: Dict[str, str], 
        nuc_A: Optional[Dict[str, str]] = None, 
        nuc_B: Optional[Dict[str, str]] = None, 
        idA: Optional[str] = None, 
        idB: Optional[str] = None, 
        con_score: str = 'y',
        ref: str = '1'
    ):
        """
        Initialize a Scoring object and calculate sequence comparison scores.
        
        Args:
            pep_A: Dictionary of aligned protein sequences for group A
            pep_B: Dictionary of aligned protein sequences for group B
            nuc_A: Dictionary of aligned nucleotide sequences for group A (optional)
            nuc_B: Dictionary of aligned nucleotide sequences for group B (optional)
            idA: ID of the reference sequence for group A (optional)
            idB: ID of the reference sequence for group B (optional)
            con_score: Whether to include conservation in score calculation ("y"/"n")
        """

        self.BLOSUM62 = Constants.BLOSUM62
        self.hydrophobicity = Constants.hydrophobicity
        self.pep_A = pep_A
        self.pep_B = pep_B
        self.nuc_A = nuc_A
        self.nuc_B = nuc_B
        self.idA = idA
        self.idB = idB
        self.con_score = con_score
        self.results = None
        self.ref = ref
        
        # Validate input sequences
        self._validate_input()
        
        # Calculate the results
        self._calculate()

    def _validate_input(self) -> None:
        """
        Validate the input sequences before processing.
        
        Checks for:
        - Empty dictionaries
        - Sequence length consistency within groups
        - Matching reference sequence IDs if provided
        
        Raises:
            ValueError: If any validation checks fail
        """
        # Check for empty dictionaries
        if not self.pep_A or not self.pep_B:
            raise ValueError("Empty sequence dictionaries provided")
        
        # Check for sequence length consistency
        pep_A_values = list(self.pep_A.values())
        pep_B_values = list(self.pep_B.values())
        
        len_A = len(pep_A_values[0])
        len_B = len(pep_B_values[0])
        
        if len_A != len_B:
            raise ValueError(f"Sequence length mismatch between groups: A ({len_A}) vs B ({len_B})")
        
        for seq in pep_A_values:
            if len(seq) != len_A:
                raise ValueError("Inconsistent sequence lengths within group A")
        
        for seq in pep_B_values:
            if len(seq) != len_B:
                raise ValueError("Inconsistent sequence lengths within group B")
        
        # Check if provided reference sequence IDs exist
        if self.idA and self.idA not in self.pep_A:
            raise ValueError(f"Reference sequence ID {self.idA} not found in group A")
        
        if self.idB and self.idB not in self.pep_B:
            raise ValueError(f"Reference sequence ID {self.idB} not found in group B")
        
        # If nucleotide sequences are provided, check if they match protein sequences
        if self.nuc_A and self.nuc_B:
            # Check if the keys match
            if set(self.pep_A.keys()) != set(self.nuc_A.keys()):
                raise ValueError("Mismatch between protein and nucleotide sequence IDs in group A")
                
            if set(self.pep_B.keys()) != set(self.nuc_B.keys()):
                raise ValueError("Mismatch between protein and nucleotide sequence IDs in group B")
            
            # Check length consistency (nucleotide should be 3x protein length if aligned)
            for key in self.pep_A:
                pep_len = len(self.pep_A[key].replace("-", ""))
                nuc_len = len(self.nuc_A[key].replace("-", ""))
                if nuc_len != pep_len * 3:
                    raise ValueError(f"Nucleotide sequence length for {key} in group A is not 3x the protein length")
            
            for key in self.pep_B:
                pep_len = len(self.pep_B[key].replace("-", ""))
                nuc_len = len(self.nuc_B[key].replace("-", ""))
                if nuc_len != pep_len * 3:
                    raise ValueError(f"Nucleotide sequence length for {key} in group B is not 3x the protein length")

    def _calculate(self) -> pd.DataFrame:
        """
        Perform the comparison calculation and store the result.
        
        Returns:
            Array containing the calculated comparison results
        """
        self.results = self._compare_seq(
            self.pep_A, self.pep_B, self.nuc_A, self.nuc_B, 
            self.idA, self.idB
        )
        return self.results

    def _define_RefSeq(self, idA: Optional[str] = None, idB: Optional[str] = None) -> Tuple[str, str]:
        """
        Establish reference sequences based on provided IDs or use the first sequence available.
        
        Args:
            idA: ID of the reference sequence for group A (optional)
            idB: ID of the reference sequence for group B (optional)
            
        Returns:
            Tuple containing reference sequences (ref_seqA, ref_seqB)
        """
        # Use provided ID if available, otherwise use the first sequence
        if idA and idA in self.pep_A:
            ref_seqA = self.pep_A[idA]
        else:
            ref_seqA = next(iter(self.pep_A.values()))
            self.idA = next(iter(self.pep_A.keys()))
            
        if idB and idB in self.pep_B:
            ref_seqB = self.pep_B[idB]
        else:
            ref_seqB = next(iter(self.pep_B.values()))
            self.idB = next(iter(self.pep_B.keys()))
            
        return ref_seqA, ref_seqB

    def _update_score(self, aa1: str, aa2: str, score: float) -> float:
        """
        Update the score based on amino acid comparison.
        
        The score is increased based on:
        - BLOSUM62 substitution matrix score
        - Hydrophobicity difference between amino acids
        - No score change for gaps or identical amino acids
        
        Args:
            aa1: Amino acid from sequence A
            aa2: Amino acid from sequence B
            score: Current score value
            
        Returns:
            Updated score value
        """
        # No score change for gaps or identical amino acids
        if aa1 == "-" or aa2 == "-" or aa1 == aa2:
            return score
        
        # Get BLOSUM62 score for the amino acid pair
        blosum_score = self.BLOSUM62.get((aa1, aa2), self.BLOSUM62.get((aa2, aa1), -5))
        
        # Calculate hydrophobicity difference
        if aa1 in self.hydrophobicity and aa2 in self.hydrophobicity:
            chem_diff = abs(self.hydrophobicity[aa1] - self.hydrophobicity[aa2])
            # Increase score based on BLOSUM62 score and hydrophobicity difference
            score += abs(blosum_score) * (1 + chem_diff)
        else:
            # Fallback if hydrophobicity values are missing
            score += abs(blosum_score)
            
        return score

    def _calculate_score(
        self,
        seq_array_A: List[str],
        seq_array_B: List[str],
        most_common_A: str,
        most_common_B: str,
        conservation_A: float,
        conservation_B: float,
    ) -> float:
        """
        Calculate score based on differences between group A and group B amino acids,
        prioritizing identical matches and then sorting remaining by hydrophobicity.

        Args:
            seq_array_A: Array of amino acids at current position for group A
            seq_array_B: Array of amino acids at current position for group B
            most_common_A: Most common amino acid in group A at this position
            most_common_B: Most common amino acid in group B at this position
            conservation_A: Conservation rate for group A at this position
            conservation_B: Conservation rate for group B at this position

        Returns:
            Position-specific score

        Raises:
            ValueError: If empty sequence arrays are provided
        """
        # Check for empty sequence arrays
        if not seq_array_A or not seq_array_B:
            # Or would you like to return 0 here? The original logic raised an error.
            raise ValueError("Empty sequence arrays")

        score = 0.0

        # 1. Count amino acids
        counter_A = collections.Counter(seq_array_A)
        counter_B = collections.Counter(seq_array_B)

        # 2. Phase: Match identical amino acids
        # We iterate over the amino acids in group A and look for identical ones in B
        # It's important to iterate over a copy of the keys when modifying the Counter object
        for aa in list(counter_A.keys()):
            if aa in counter_B:
                # How many identical pairs can we form?
                num_identical_pairs = min(counter_A[aa], counter_B[aa])

                # Match identical pairs and update score
                for _ in range(num_identical_pairs):
                    score = self._update_score(aa, aa, score)

                # Remove the matched amino acids from the Counters
                counter_A[aa] -= num_identical_pairs
                counter_B[aa] -= num_identical_pairs

        # 3. Phase: Prepare and sort remaining amino acids by hydrophobicity
        remaining_A = list(counter_A.elements()) # Convert Counter back to a list of remaining AAs
        remaining_B = list(counter_B.elements())

        # Sort by hydrophobicity (descending)
        # We use .get(aa, 0) in case an AA is not in the dictionary for any reason
        sorted_remaining_A = sorted(
            remaining_A,
            key=lambda aa: Constants.hydrophobicity.get(aa, 0),
            reverse=True
        )
        sorted_remaining_B = sorted(
            remaining_B,
            key=lambda aa: Constants.hydrophobicity.get(aa, 0),
            reverse=True
        )

        # 4. Phase: Match the remaining, hydrophobicity-sorted amino acids
        len_remaining_A = len(sorted_remaining_A)
        len_remaining_B = len(sorted_remaining_B)
        max_len_remaining = max(len_remaining_A, len_remaining_B)

        for j in range(max_len_remaining):
            aa1, aa2 = None, None

            if j < len_remaining_A and j < len_remaining_B:
                # Both lists still have elements
                aa1 = sorted_remaining_A[j]
                aa2 = sorted_remaining_B[j]
            elif j < len_remaining_A:
                # Group B has no remaining elements, compare remaining A with most_common_B
                aa1 = sorted_remaining_A[j]
                aa2 = most_common_B
            elif j < len_remaining_B:
                 # Group A has no remaining elements, compare remaining B with most_common_A
                aa1 = most_common_A
                aa2 = sorted_remaining_B[j]

            # Only update if aa1 and aa2 have been set (should always be the case in the loop)
            if aa1 is not None and aa2 is not None:
                 score = self._update_score(aa1, aa2, score)


        # 5. Adjust score by conservation rates (as in the original)
        if hasattr(self, 'con_score') and self.con_score.upper() in ["Y", "YES"]:
             # Ensure that conservation_A and conservation_B are not None or 0,
             # in case this could lead to unexpected results.
             if conservation_A is not None and conservation_B is not None:
                 score = score * conservation_A * conservation_B
             else:
                 # Optional warning or error handling if Conservation is None
                 print("Warning: Conservation rates are None. Skipping conservation weighting.")


        return score

    def _update_Common_Conservation(
        self, 
        array: List[str],
        return_conservation: bool = True,
        return_frequencies: bool = True
        ) -> Tuple[str, Optional[float], Optional[Dict[str, float]]]:
        """
        Update the most common element list and optionally the conservation rate list
        and frequency dictionaries list.
        
        Args:
            array: List of elements at the current position
            most_common_list: List to append the most common element to
            conservation_list: List to append the conservation rate to (optional)
            frequency_dicts_list: List to append the frequency dictionary to (optional)
            return_conservation: Whether to calculate and return conservation rates
            return_frequencies: Whether to calculate and return frequency dictionaries
            
        Returns:
            Depending on parameters:
            - most_common_list only
            - (most_common_list, conservation_list)
            - (most_common_list, frequency_dicts_list)
            - (most_common_list, conservation_list, frequency_dicts_list)
        """
        # Count occurrences of each element
        counter = collections.Counter(array)
        max_count = max(counter.values()) if counter else 0
        
        # Find all elements with the maximum count (in case of ties)
        most_common_elements = [key for key, count in counter.items() if count == max_count]
        
        # Randomly choose one of the most common elements
        most_common = random.choice(most_common_elements) if most_common_elements else ""
        #most_common_list.append(most_common)
        
        # Calculate conservation rate if requested
        if return_conservation:
            # if conservation_list is None:
            #     conservation_list = []
                
            if len(array) == 0:
                conservation_rate = 0.0
            else:
                conservation_rate = float(max_count / len(array))
        
        # Calculate frequency dictionary if requested
        if return_frequencies:
            # if frequency_dicts_list is None:
            #     frequency_dicts_list = []
                
            # Create frequency dictionary with relative frequencies
            freq_dict = {}
            total = len(array) if len(array) > 0 else 1  # Avoid division by zero
            
            for element, count in counter.items():
                freq_dict[element] = float(count / total)
                
            # frequency_dicts_list.append(freq_dict)

        return most_common, conservation_rate, freq_dict
    
    def _remove_gaps(self, data_dict: Dict[str, Any]) -> Dict[str, Any]:
        """
        Remove positions where both groups have gaps in the most common amino acids.
        
        Args:
            data_dict: Dictionary containing position-specific data
            
        Returns:
            Dictionary with gap positions removed
        """
        # Create a copy to avoid modifying the original dictionary
        result_dict = data_dict.copy()
        
        # Check if nucleotide data is present
        has_triplets = 'MostCommonTripletA' in result_dict and 'MostCommonTripletB' in result_dict
        
        seq_length = len(result_dict['ref_seqA'])
        
        # Lists to track which positions to keep
        positions_to_keep = []

        # First pass: identify positions to keep
        for k in range(seq_length):
            if result_dict['ref_seqA'][k] != "-" or result_dict['ref_seqB'][k] != "-":
                positions_to_keep.append(k)
        
        # Create new arrays with only the positions to keep
        for key in ['MostCommonA', 'MostCommonB', 'ConservationA', 'ConservationB', 'FrequenciesA', 'FrequenciesB', 'Score']:
            result_dict[key] = [result_dict[key][i] for i in positions_to_keep]
        
        # Handle nucleotide data if present
        if has_triplets:
            result_dict['MostCommonTripletA'] = [result_dict['MostCommonTripletA'][i] for i in positions_to_keep]
            result_dict['MostCommonTripletB'] = [result_dict['MostCommonTripletB'][i] for i in positions_to_keep]
            result_dict['PointMutations'] = [result_dict['PointMutations'][i] for i in positions_to_keep]
        
        # Create new reference sequences
        result_dict['ref_seqA'] = ''.join(result_dict['ref_seqA'][i] for i in positions_to_keep)
        result_dict['ref_seqB'] = ''.join(result_dict['ref_seqB'][i] for i in positions_to_keep)
        
        return result_dict

    def _define_reference_positions(self, ref_seq: str) -> List[Union[int, str]]:
        """
        Create a list of reference positions, handling gaps in the reference sequence.
        
        Args:
            ref_seq: Reference sequence string
            
        Returns:
            List of reference positions (integers for amino acids, "nan" for gaps)
        """
        reference_positions = []
        pos = 0
        
        for aa in ref_seq:
            if aa == "-":
                reference_positions.append("nan")
            else:
                pos += 1
                reference_positions.append(pos)
                
        return reference_positions
    
    def _triplets_with_point_mutation(self, first_triplets: Tuple[str], second_triplets: Tuple[str]) -> List[str]:
        """
        Takes the most_common_triplets arrays and searches for single character differences.
        
        Args:
            first_triplets: Tuple containing the most common triplets in the first group
            second_triplets: Tuple containing the most common triplets in the second group
            
        Returns:
            List with information whether each aa difference could be due to point mutation
        """
        point_mutations: List[str] = []
        for i in range(len(first_triplets)):
            differences = 0
            for j in range(len(first_triplets[i])):
                if first_triplets[i][j] != second_triplets[i][j]:
                    differences += 1
            if differences == 1:
                point_mutations.append('Y')
            else:
                point_mutations.append('N')
        return point_mutations

    def _compare_seq(
        self, 
        pep_A: Dict[str, str], 
        pep_B: Dict[str, str], 
        nuc_A: Optional[Dict[str, str]] = None, 
        nuc_B: Optional[Dict[str, str]] = None, 
        idA: Optional[str] = None, 
        idB: Optional[str] = None, 
    ) -> pd.DataFrame:
        """
        Compare two groups of sequences and calculate position-specific scores.
        
        This method:
        1. Analyzes each position in the aligned sequences
        2. Determines the most common amino acid at each position for both groups
        3. Calculates conservation rates
        4. Computes similarity scores
        5. Removes positions where both groups have gaps
        6. Normalizes scores
        
        Args:
            pep_A: Dictionary of aligned protein sequences for group A
            pep_B: Dictionary of aligned protein sequences for group B
            nuc_A: Dictionary of aligned nucleotide sequences for group A (optional)
            nuc_B: Dictionary of aligned nucleotide sequences for group B (optional)
            idA: ID of the reference sequence for group A (optional)
            idB: ID of the reference sequence for group B (optional)
            con_score: Whether to include conservation in score calculation ("y"/"n")
            
        Returns:
            Array containing the calculated comparison results
        """
        # Initialize result data structures
        most_common_A = []
        most_common_B = []
        frequencies_A: List[Dict[str, float]] = []
        frequencies_B: List[Dict[str, float]] = []
        conservation_A = []
        conservation_B = []
        most_common_triplet_A = []
        most_common_triplet_B = []
        point_mutations = []
        scores = []

        reference_positions_A = []
        reference_positions_B = []
        
        # Reference sequence length (assumes all sequences are aligned and have same length)
        seq_length = len(next(iter(pep_A.values())))
        
        # Process each position in the aligned sequences
        nuc_pos = 0
        for pos in range(seq_length):
            # Extract amino acids at current position
            aa_array_A = [seq[pos] for seq in pep_A.values()]
            aa_array_B = [seq[pos] for seq in pep_B.values()]
            
            # Update most common amino acids and conservation rates
            # most_common_A, conservation_A, frequencies_A = self._update_Common_Conservation(
            #     aa_array_A, most_common_A, conservation_A, frequencies_A
            # )
            most_common, conservation, frequencies = self._update_Common_Conservation(aa_array_A)
            most_common_A.append(most_common)
            conservation_A.append(conservation)
            frequencies_A.append(frequencies)
            
            # most_common_B, conservation_B, frequencies_B = self._update_Common_Conservation(
            #     aa_array_B, most_common_B, conservation_B, frequencies_B
            # )
            most_common, conservation, frequencies = self._update_Common_Conservation(aa_array_B)
            most_common_B.append(most_common)
            conservation_B.append(conservation)
            frequencies_B.append(frequencies)
            
            # Process nucleotide sequences if provided
            if nuc_A and nuc_B:
                # Extract nucleotide triplets at current position
                triplet_array_A = [seq[nuc_pos:nuc_pos+3] for seq in nuc_A.values()]
                triplet_array_B = [seq[nuc_pos:nuc_pos+3] for seq in nuc_B.values()]
                
                # Update most common triplets
                # most_common_triplet_A = self._update_Common_Conservation(
                #     triplet_array_A, most_common_triplet_A, return_conservation=False, return_frequencies=False
                # )
                most_common = self._update_Common_Conservation(triplet_array_A)
                most_common_triplet_A.append(most_common)

                
                # most_common_triplet_B = self._update_Common_Conservation(
                #     triplet_array_B, most_common_triplet_B, return_conservation=False, return_frequencies=False
                # )
                most_common = self._update_Common_Conservation(triplet_array_B)
                most_common_triplet_B.append(most_common)
                
                # Update nucleotide position (3 nucleotides per amino acid)
                nuc_pos += 3
            
            # Calculate position-specific score
            score = self._calculate_score(
                aa_array_A, aa_array_B, 
                most_common_A[pos], most_common_B[pos],
                conservation_A[pos], conservation_B[pos]
            )
            
            scores.append(score)
        
        # Get reference sequences
        ref_seq_A, ref_seq_B = self._define_RefSeq(idA, idB)

        # Get point mutations
        point_mutations = self._triplets_with_point_mutation(most_common_triplet_A, most_common_triplet_B)

        # Create data dictionary for gap removal
        data_dict = {
            'MostCommonA': most_common_A,
            'MostCommonB': most_common_B,
            'ConservationA': conservation_A,
            'ConservationB': conservation_B,
            'FrequenciesA': frequencies_A,
            'FrequenciesB': frequencies_B,
            'Score': scores,
            'ref_seqA': ref_seq_A,
            'ref_seqB': ref_seq_B
        }

        # Add nucleotide data if available
        if nuc_A and nuc_B:
            data_dict['MostCommonTripletA'] = most_common_triplet_A
            data_dict['MostCommonTripletB'] = most_common_triplet_B
            data_dict['PointMutations'] = point_mutations
        
        # Remove positions where both groups have gaps
        cleaned_data = self._remove_gaps(data_dict)
        
        # Extract cleaned data
        most_common_A = cleaned_data['MostCommonA']
        most_common_B = cleaned_data['MostCommonB']
        conservation_A = cleaned_data['ConservationA']
        conservation_B = cleaned_data['ConservationB']
        frequencies_A = cleaned_data['FrequenciesA']
        frequencies_B = cleaned_data['FrequenciesB']
        scores = cleaned_data['Score']
        ref_seq_A = cleaned_data['ref_seqA']
        ref_seq_B = cleaned_data['ref_seqB']

        # Create position indices and reference positions
        positions = list(range(1, len(scores) + 1))
        reference_positions_A = self._define_reference_positions(ref_seq_A)
        reference_positions_B = self._define_reference_positions(ref_seq_B)
        
        # Normalize scores
        if scores:
            lenA = len(self.pep_A.items())
            lenB = len(self.pep_B.items())
            if lenA > lenB:
                normalized_scores = [score / lenA for score in scores]
            else:
                normalized_scores = [score / lenB for score in scores]
            # score_max = max(scores) if max(scores) > 0 else 1.0
            # normalized_scores = [score / score_max for score in scores]
        else:
            normalized_scores = []
        
        # return result_array
        if nuc_A and nuc_B:
            result_df = pd.DataFrame({
                "POSITION": positions,
                "REF_POS_A": reference_positions_A,
                "REF_POS_B": reference_positions_B,
                "SCORE": normalized_scores,
                "COMMON_A": most_common_A,
                "CONSERVATION_A": conservation_A,
                "FREQUENCIES_A": frequencies_A,
                "COMMON_B": most_common_B,
                "CONSERVATION_B": conservation_B,
                "FREQUENCIES_B": frequencies_B,
                "TRIPLET_A": cleaned_data["MostCommonTripletA"],
                "TRIPLET_B": cleaned_data["MostCommonTripletB"],
                "POINT_MUTATIONS": cleaned_data["PointMutations"]
            })
        else:
            result_df = pd.DataFrame({
                "POSITION": positions,
                "REF_POS_A": reference_positions_A,
                "REF_POS_B": reference_positions_B,
                "SCORE": normalized_scores,
                "COMMON_A": most_common_A,
                "CONSERVATION_A": conservation_A,
                "FREQUENCIES_A": frequencies_A,
                "COMMON_B": most_common_B,
                "CONSERVATION_B": conservation_B,
                "FREQUENCIES_B": frequencies_B
            })

        return result_df

# endregion

# region Output
class Output:
    """
    A class for processing, saving, and visualizing comparison results between two sequence groups.
    
    This class handles the processing of score dataframes generated from sequence comparisons,
    provides visualization capabilities using matplotlib or plotly, and allows for
    saving results to CSV files.
    """
    
    def __init__(
        self,
        result_df: pd.DataFrame,
        show: str,
        top: int = 20, 
        lb: Optional[int] = None, 
        ub: Optional[int] = None,
        seqpos: Optional[Tuple[int]] = None,
        group_a: str = "A", 
        group_b: str = "B", 
        ref: str = "1",
        idA: str = "First Sequence of Group 1",
        idB: str = "First Sequence of Group 2",
        directory: str = None
    ) -> None:
        """
        Initialize the Output object with result data and display parameters.
        
        Args:
            result_df: DataFrame containing position, reference positions, common residues, 
                      conservation values and scores.
            top: Number of top-scoring positions to show when not using bounds.
            lb: Lower bound for filtering positions (based on reference sequence).
            ub: Upper bound for filtering positions (based on reference sequence).
            group_a: Label for the first sequence group.
            group_b: Label for the second sequence group.
            ref: Which reference sequence to use ("1" or "2").
            idA: Reference-ID from first group.
            idB: Reference-ID from second group.
        """
        self.result_df = result_df
        self.show = show
        self.top = top
        self.lb = lb
        self.ub = ub
        self.seqpos = seqpos
        self.group_a = group_a
        self.group_b = group_b
        self.ref = ref
        self.idA = idA
        self.idB = idB
        self.directory = directory

        # Column names for better readability
        self.POSITION = 'POSITION'
        self.REF_POS_A = 'REF_POS_A'
        self.REF_POS_B = 'REF_POS_B'
        self.SCORE = 'SCORE'
        self.COMMON_A = 'COMMON_A'
        self.CONSERVATION_A = 'CONSERVATION_A'
        self.FREQUENCIES_A = 'FREQUENCIES_A'
        self.COMMON_B = 'COMMON_B'
        self.CONSERVATION_B = 'CONSERVATION_B'
        self.FREQUENCIES_B = 'FREQUENCIES_B'
        self.TRIPLET_A = 'TRIPLET_A'
        self.TRIPLET_B = 'TRIPLET_B'
        self.POINT_MUTATIONS = 'POINT_MUTATIONS'
        
    def save(self) -> None:
            #def save(self, directory: Optional[str] = None) -> None:
        """
        Saves the result dataframe as a CSV file in the specified directory.
        
        Uses the parameters specified during initialization to save results
        to a 'CoSeD_output.csv' file in the specified directory.
        
        Args:
            directory: Directory to save the CSV file. If None, uses the directory
                      specified during initialization.
        
        Raises:
            ValueError: If no directory specified for saving output.
        """
        if not self.directory:
            raise ValueError("No directory specified for saving output")
            
        # Ensure directory exists
        os.makedirs(self.directory, exist_ok=True)
            
        # Create a copy of the dataframe for saving
        result_df_csv = self.result_df.copy()

        # Create the header with appropriate column names
        if self.TRIPLET_A in self.result_df.columns:
            columns = {
                self.POSITION: f"Position",
                self.REF_POS_A: f"Positions_in_{self.idA}",
                self.REF_POS_B: f"Positions_in_{self.idB}",
                self.SCORE: f"Score",
                self.COMMON_A: f"Most_Common_{self.group_a}",
                self.CONSERVATION_A: f"Conservation_{self.group_a}",
                self.COMMON_B: f"Most_Common_{self.group_b}",
                self.CONSERVATION_B: f"Conservation_{self.group_b}",
                self.TRIPLET_A: f"Most_Common_Triplet_{self.group_a}",
                self.TRIPLET_B: f"Most_Common_Triplet_{self.group_b}",
                self.POINT_MUTATIONS: f"Point_Mutations",
                self.FREQUENCIES_A: f"AA_frequencies_{self.group_a}",
                self.FREQUENCIES_B: f"AA_frequencies_{self.group_b}"
            }
        else:
            columns = {
                self.POSITION: f"Position",
                self.REF_POS_A: f"Positions_in_{self.idA}_(first_group)",
                self.REF_POS_B: f"Positions_in_{self.idB}",
                self.SCORE: f"Score",
                self.COMMON_A: f"Most_Common_{self.group_a}",
                self.CONSERVATION_A: f"Conservation_{self.group_a}",
                self.COMMON_B: f"Most_Common_{self.group_b}",
                self.CONSERVATION_B: f"Conservation_{self.group_b}",
                self.FREQUENCIES_A: f"AA_frequencies_{self.group_a}",
                self.FREQUENCIES_B: f"AA_frequencies_{self.group_b}"
            }
        
        # Rename columns for saving
        result_df_csv = result_df_csv.rename(columns=columns)
        
        # Save to CSV file
        output_path = f"{self.directory}/CoSeD_output.csv"
        result_df_csv.to_csv(output_path, index=False)

        print(f"CSV saved to {output_path}")

    @staticmethod
    def check_display_tools() -> Dict[str, bool]:
        """
        Checks whether required visualization packages are installed.
        
        Returns:
            Dictionary mapping tool names to booleans indicating availability.
        """
        tools = ["matplotlib", "plotly"]
        installed_tools = {tool: importlib.util.find_spec(tool) is not None for tool in tools}
        return installed_tools
    
    def _round_numerical_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Helper method to round numerical columns to two decimal places.
        
        Args:
            df: The pandas DataFrame with columns to round.
            
        Returns:
            DataFrame with rounded numerical columns.
        """
        # Create a copy to avoid modifying the original
        result = df.copy()
        
        # Round numerical columns
        result[self.SCORE] = result[self.SCORE].round(2)
        result[self.CONSERVATION_A] = result[self.CONSERVATION_A].round(2)
        result[self.CONSERVATION_B] = result[self.CONSERVATION_B].round(2)
        
        return result

    def _generate_display_df(self) -> pd.DataFrame:
        """
        Generates a filtered and processed dataframe for visualization.
        
        Applies filtering based on specific sequence positions, position ranges,
        or top scores, and performs formatting operations to prepare data for display.
        
        Returns:
            Processed pandas DataFrame ready for visualization.
        """
        result_df = self.result_df
        
        # Apply filtering based on specific sequence positions if provided
        if len(self.seqpos) > 0:
            # Column for reference sequence (REF_POS_B for ref="2", REF_POS_A for ref="1") 
            ref_col = self.REF_POS_B if self.ref == "2" else self.REF_POS_A
            
            # Converting ref_col to numeric values
            result_df[ref_col] = pd.to_numeric(result_df[ref_col], errors='coerce')
            
            # Filter based on specific sequence positions from the tuple
            display_df = result_df[result_df[ref_col].isin(self.seqpos)].copy()
            
            # Process reference position column
            # Convert to numeric, replace non-numeric values with NaN
            display_df[self.REF_POS_A] = pd.to_numeric(display_df[self.REF_POS_A], errors='coerce')
            
            # Remove rows with NaN in reference positions
            display_df = display_df.dropna(subset=[self.REF_POS_A])
            
            # Format reference positions: convert to int and then to string
            display_df[self.REF_POS_A] = display_df[self.REF_POS_A].astype(int).astype(str)
        
        # Apply position-based filtering if bounds are provided
        elif self.lb is not None and self.ub is not None:
            lb_val = max(1.0, float(self.lb))
            
            # Get maximum position from the result dataframe
            pos_max = result_df[self.POSITION].max()
            ub_val = min(float(self.ub), pos_max)
            
            # Column for reference sequence (REF_POS_B for ref="2", REF_POS_A for ref="1") 
            ref_col = self.REF_POS_B if self.ref == "2" else self.REF_POS_A
            
            # Converting ref_col to numeric values
            result_df[ref_col] = pd.to_numeric(result_df[ref_col], errors='coerce')
            
            # Filter based on reference positions
            display_df = result_df[(result_df[ref_col] >= lb_val) & (result_df[ref_col] <= ub_val)].copy()
            
            # Process reference position column
            # Convert to numeric, replace non-numeric values with NaN
            display_df[self.REF_POS_A] = pd.to_numeric(display_df[self.REF_POS_A], errors='coerce')
            
            # Remove rows with NaN in reference positions
            display_df = display_df.dropna(subset=[self.REF_POS_A])
            
            # Format reference positions: convert to int and then to string
            display_df[self.REF_POS_A] = display_df[self.REF_POS_A].astype(int).astype(str)

        else:
            # First select the records with the highest scores
            top_scores_df = result_df.nlargest(self.top, self.SCORE).copy()
            
            # Then sort by reference position
            ref_col = self.REF_POS_B if self.ref == "2" else self.REF_POS_A
            
            # Converting the reference position to numeric values ​​for sorting
            top_scores_df[ref_col] = pd.to_numeric(top_scores_df[ref_col], errors='coerce')
            
            # Sort by reference position
            display_df = top_scores_df.sort_values(by=ref_col).copy()
        
        # Round numerical values to two decimals for better display
        display_df = self._round_numerical_columns(display_df)
        
        return display_df

    def _generate_point_mutation_df(self, top: bool = True) -> pd.DataFrame:
        """
        Generates a filtered and processed dataframe for visualization of point mutations.
        
        Applies filtering based on positions of point mutations and top scores, 
        and performs formatting operations to prepare data for display.
        
        Args:
            top: If True, only return the top-scoring point mutations based on self.top.
        
        Returns:
            Processed pandas DataFrame ready for visualization.
        """
        result_df = self.result_df
        
        # Check if point mutation column exists
        if self.POINT_MUTATIONS not in result_df.columns:
            raise ValueError("Point mutation data not available in the result dataframe")
            
        # Choose lines with point mutations and different amino acids
        display_df = result_df[
            (result_df[self.POINT_MUTATIONS] == 'Y') & 
            (result_df[self.COMMON_A] != result_df[self.COMMON_B])
        ].copy()
        
        # Sort by score in descending order and take top rows if requested
        if top:
            display_df = display_df.sort_values(by=self.SCORE, ascending=False).head(self.top)
        
        # Round numerical values
        display_df = self._round_numerical_columns(display_df)
        
        return display_df

    def _get_plot_data(self, display_df: pd.DataFrame) -> Tuple[pd.Series, ...]:
        """
        Extract data from the display dataframe needed for plotting.
        
        Args:
            display_df: Processed dataframe for display.
            
        Returns:
            Tuple containing extracted data: (x, y_up, y_down, labels_a, labels_b, 
                                            score, has_triplets, triplet_a, triplet_b, x_labels)
        """
        # Determine x-axis based on the chosen reference sequence
        ref_col = self.REF_POS_B if self.ref == "2" else self.REF_POS_A
        
        # Store original reference positions for labels
        x_labels = display_df[ref_col]
        
        # Create consecutive indices for bars when using 'top' mode
        if self.lb is None or self.ub is None:
            # We're in 'top' mode - use consecutive indices (0, 1, 2...)
            x = pd.Series(range(len(display_df)))
        else:
            # We're in 'bounds' mode - use actual reference positions
            x = display_df[ref_col]
        
        # Extract data for plotting
        y_up = display_df[self.CONSERVATION_A].abs()
        y_down = -display_df[self.CONSERVATION_B].abs()

        labels_a = display_df[self.COMMON_A]
        labels_b = display_df[self.COMMON_B]
        score = display_df[self.SCORE]
        
        # Check if triplet data is available
        has_triplets = self.TRIPLET_A in display_df.columns
        if has_triplets:
            triplet_a = display_df[self.TRIPLET_A]
            triplet_b = display_df[self.TRIPLET_B]
        else:
            triplet_a = pd.Series([""] * len(labels_a))
            triplet_b = pd.Series([""] * len(labels_b))
        
        return x, y_up, y_down, labels_a, labels_b, score, has_triplets, triplet_a, triplet_b, x_labels

    def _display_matplotlib(self, display_df: Optional[pd.DataFrame] = None, pm_plot: str = 'Y', id: str = "1") -> None:
        """
        Displays a bar plot of the result dataframe using matplotlib.
        
        Creates a visualization showing conservation values for both groups, along
        with annotations for common residues and scores.
        
        Args:
            display_df: Pre-processed dataframe for display. If None, generated automatically.
            pm_plot: Whether to create a point mutation plot. If None, uses self.point_mutations.
            id: ID of reference sequence.
        """
        if plt is None:
            raise ImportError("Matplotlib is not installed. Please install matplotlib to use this feature.")
            
        # Use provided display dataframe or generate one
        if display_df is None:
            display_df = self._generate_display_df()
            
        # Extract plot data with additional x_labels
        x, y_up, y_down, labels_a, labels_b, score, has_triplets, triplet_a, triplet_b, x_labels = self._get_plot_data(display_df)
        
        # Create the bar plot
        fig, ax = plt.subplots(figsize=(13, 8))
        ax.bar(x, y_up, color='darkorange', label=self.group_a)
        ax.bar(x, y_down, color='darkgreen', label=self.group_b)
        
        # Annotate bars
        for i in range(len(x)):
            annotation_a = f"{labels_a.iloc[i]}"
            annotation_b = f"{labels_b.iloc[i]}"
            
            if has_triplets:
                annotation_a += f"\n\n{triplet_a.iloc[i]}"
                annotation_b += f"\n\n{triplet_b.iloc[i]}"
                
            ax.text(x.iloc[i], y_up.iloc[i] - 0.05, annotation_a, ha='center', va='top', fontsize=12, color='white')
            ax.text(x.iloc[i], y_down.iloc[i] + 0.05, annotation_b, ha='center', va='bottom', fontsize=12, color='white')
            ax.text(x.iloc[i], 0.05, f"Score\n{score.iloc[i]}", ha='center', va='bottom', fontsize=8, color='white')
        
        # Configure plot layout
        ax.set_xlabel(f'Indices from {id}')
        ax.set_ylabel('Conservation')
        ax.set_title('Plot of AA Conservation')
        
        # Set x-ticks to show original reference positions when in 'top' mode
        if self.lb is None or self.ub is None:
            ax.set_xticks(x)
            ax.set_xticklabels(x_labels, rotation=45, ha='right')
        
        ax.axhline(0, color='black', linewidth=1)  # Zero line
        ax.set_ylim(-1, 1)
        ax.legend()
        plt.subplots_adjust(bottom=0.2)

        # Create point mutation plot if requested and data available
        if pm_plot.upper() in ['Y', 'YES'] and has_triplets:
            try:
                display_df_pm = self._generate_point_mutation_df(top=True)
                if len(display_df_pm) == 0:
                    print("No point mutations found.")
                    return
                    
                # Extract plot data for point mutations with additional x_labels
                x, y_up, y_down, labels_a, labels_b, score, _, triplet_a, triplet_b, x_labels = self._get_plot_data(display_df_pm)

                # Create the bar plot
                fig, ax = plt.subplots(figsize=(13, 8))
                ax.bar(x, y_up, color='darkcyan', label=self.group_a)
                ax.bar(x, y_down, color='pink', label=self.group_b)
                
                # Annotate bars
                for i in range(len(x)):
                    annotation_a = f"{labels_a.iloc[i]}\n\n{triplet_a.iloc[i]}"
                    annotation_b = f"{labels_b.iloc[i]}\n\n{triplet_b.iloc[i]}"
                    ax.text(x.iloc[i], y_up.iloc[i] - 0.05, annotation_a, ha='center', va='top', fontsize=12, color='white')
                    ax.text(x.iloc[i], y_down.iloc[i] + 0.05, annotation_b, ha='center', va='bottom', fontsize=12, color='white')
                    ax.text(x.iloc[i], 0.05, f"Score\n{score.iloc[i]}", ha='center', va='bottom', fontsize=8, color='white')
                
                # Configure plot layout
                ax.set_xlabel(f'Indices from {id}')
                ax.set_ylabel('Conservation')
                ax.set_title('Plot of point mutations')
                
                # Set x-ticks to show original reference positions when in 'top' mode
                if self.lb is None or self.ub is None:
                    ax.set_xticks(x)
                    ax.set_xticklabels(x_labels, rotation=45, ha='right')
                
                ax.axhline(0, color='black', linewidth=1)  # Zero line
                ax.set_ylim(-1, 1)
                ax.legend()
                plt.subplots_adjust(bottom=0.2)
                
            except ValueError as e:
                print(f"Could not generate point mutation plot: {e}")

        # Save Matplotlib figure
        save_path = os.path.join(self.directory, "matplotlib_plot.png")
        plt.savefig(save_path, dpi=300, bbox_inches='tight')    

        if self.show.upper() in ["YES" or "Y"]:
            plt.show() # Show the plot

    def _display_plotly(self, display_df: Optional[pd.DataFrame] = None, nuc_plot: str = 'Y', pm_plot: str = 'Y', id: str = 'Reference Sequence') -> None:
        """
        Displays interactive bar plots of the result dataframe using Plotly.
        
        Creates interactive visualizations showing conservation values for both groups,
        with hover information for common residues and scores.
        
        Args:
            display_df: Pre-processed dataframe for display. If None, generated automatically.
            nuc_plot: Whether to create a separate nucleotide plot.
            pm_plot: Whether to create a point mutation plot. If None, uses self.point_mutations.
            id: ID of reference sequence.
        """
        if go is None:
            raise ImportError("Plotly is not installed. Please install plotly to use this feature.")
            
        # Use provided display dataframe or generate one
        if display_df is None:
            display_df = self._generate_display_df()
            
        # Use provided pm_plot preference or fall back to initialization value
        show_pm_plot = pm_plot
            
        # Extract plot data with additional x_labels
        x, y_up, y_down, labels_a, labels_b, score, has_triplets, triplet_a, triplet_b, x_labels = self._get_plot_data(display_df)

        # Create amino acid conservation plot
        self._create_plotly_figure(
            x=x, 
            y_up=y_up, 
            y_down=y_down, 
            text_a=labels_a, 
            text_b=labels_b, 
            hover_a=triplet_a if has_triplets else labels_a, 
            hover_b=triplet_b if has_triplets else labels_b, 
            score=score, 
            color_a='darkorange', 
            color_b='darkgreen', 
            title='Plot of AA Conservation',
            plotly_id=id,
            x_labels=x_labels
        )
        
        # Create separate triplet conservation plot if requested and data available
        if nuc_plot.upper() in ['Y', 'YES'] and has_triplets:
            self._create_plotly_figure(
                x=x, 
                y_up=y_up, 
                y_down=y_down, 
                text_a=triplet_a, 
                text_b=triplet_b, 
                hover_a=labels_a, 
                hover_b=labels_b, 
                score=score, 
                color_a='darkred', 
                color_b='lightblue', 
                title='Plot of Triplet Conservation',
                plotly_id=id,
                x_labels=x_labels
            )

        # Create separate point mutation plot if requested and data available
        if show_pm_plot.upper() in ['Y', 'YES'] and has_triplets:
            try:
                display_df_pm = self._generate_point_mutation_df(top=True)
                if len(display_df_pm) == 0:
                    print("No point mutations found.")
                    return
                    
                # Extract plot data for point mutations
                x, y_up, y_down, labels_a, labels_b, score, _, triplet_a, triplet_b, x_labels = self._get_plot_data(display_df_pm)

                self._create_plotly_figure(
                    x=x, 
                    y_up=y_up, 
                    y_down=y_down, 
                    text_a=labels_a, 
                    text_b=labels_b, 
                    hover_a=triplet_a, 
                    hover_b=triplet_b, 
                    score=score, 
                    color_a='darkcyan', 
                    color_b='pink', 
                    title='Plot of point mutations',
                    plotly_id=id,
                    x_labels=x_labels
                )
            except ValueError as e:
                print(f"Could not generate point mutation plot: {e}")
            
    def _create_plotly_figure(
        self, 
        x: pd.Series, 
        y_up: pd.Series, 
        y_down: pd.Series, 
        text_a: pd.Series, 
        text_b: pd.Series, 
        hover_a: pd.Series, 
        hover_b: pd.Series, 
        score: pd.Series, 
        color_a: str, 
        color_b: str, 
        title: str,
        plotly_id: str,
        x_labels: Optional[pd.Series] = None
    ) -> None:
        """
        Helper method to create and display a Plotly figure with the provided data.
        
        Args:
            x: X-axis values (positions).
            y_up: Y values for the upper bars (group A).
            y_down: Y values for the lower bars (group B).
            text_a: Text labels for group A bars.
            text_b: Text labels for group B bars.
            hover_a: Hover text for group A bars.
            hover_b: Hover text for group B bars.
            score: Score values to display.
            color_a: Color for group A bars.
            color_b: Color for group B bars.
            title: Plot title.
            plotly_id: ID of reference sequence.
            x_labels: Original reference positions for x-axis labels.
        """
        fig = go.Figure()
        
        #Add bars for group A (positive conservation)
        fig.add_trace(go.Bar(
            x=x,
            y=y_up,
            name=self.group_a,
            marker_color=color_a,
            text=text_a,
            hovertext=hover_a,
            hoverinfo='text',
            textposition='inside',
        ))
        
        # Add bars for group B (negative conservation)
        fig.add_trace(go.Bar(
            x=x,
            y=y_down,
            name=self.group_b,
            marker_color=color_b,
            text=text_b,
            hovertext=hover_b,
            hoverinfo='text',
            textposition='inside',
        ))
        
        # Add bars for the score (transparent with text only)
        fig.add_trace(go.Bar(
            x=x,
            y=y_up * 0 + 0.1,  # Small offset for visibility
            marker_color='rgba(0,0,0,0)',  # Transparent bars
            marker_line=dict(color='white', width=0.0),
            text=score,
            textposition='outside',
            showlegend=False,
            hoverinfo='none'
        ))
        
        # Add zero horizontal line
        fig.add_trace(go.Scatter(
            x=[x.min(), x.max()],
            y=[0, 0],
            mode='lines',
            line=dict(color='black', width=1),
            showlegend=False,
            hoverinfo='none'
        ))
        
        # Customize the layout
        fig.update_layout(
            title=title,
            xaxis_title=f'Indices from {plotly_id}',
            yaxis_title='Conservation',
            barmode='relative',
            showlegend=True,
            plot_bgcolor='white',
            height=600,
            annotations=[
                dict(
                    x=0.5,
                    y=max(y_up.max(), abs(y_down.min())) + 0.4,
                    xref='paper',
                    yref='y',
                    text='Scores calculated from BLOSUM62 and Hydrophobicity (Kyle and Doolittle scale)',
                    showarrow=False,
                    font=dict(size=14)
        )])
        
        # Set x-ticks to show original reference positions when in 'top' mode
        if self.lb is None or self.ub is None and x_labels is not None:
            fig.update_layout(
                xaxis=dict(
                    tickmode='array',
                    tickvals=list(x),
                    ticktext=list(x_labels),
                    tickangle=45
                )
            )
        # Save Plotly figure
        if importlib.util.find_spec("kaleido"):
            save_path = os.path.join(self.directory, "plotly_plot.png")
            fig.write_image(save_path, width=800, height=600, scale=3)
            print(f"Plotly Plot saved as PNG to: {save_path}")
        else:
        # except ImportError:
            html_path = os.path.join(self.directory, "plotly_plot.html")
            fig.write_html(html_path)
            print(f"Kaleido not installed. Plot saved as HTML to: {html_path}")

        if self.show.upper() in ["YES", "Y"]:
            fig.show() # Show the plot


    def _create_logo_plot(self, display_df: pd.DataFrame, use_plotly: str = 'Y') -> None:
        """
        Creates a sequence logo plot based on the frequency data in both sequence groups.
        
        The logo plot displays the frequency distribution of amino acids at each position,
        with group A displayed above the x-axis and group B below it.
        
        Args:
            use_plotly: Whether to use Plotly ('Y'/'YES') or Matplotlib for visualization.
            display_df: Pre-processed dataframe for display. If None, generated automatically.
        """
        # Check which display tools are available
        tools = self.check_display_tools()
        
        # Use provided display dataframe or generate one
        if display_df is None:
            display_df = self._generate_display_df()
        
        # Extract plot data with additional x_labels
        x, _, _, _, _, score, _, _, _, x_labels = self._get_plot_data(display_df)
        
        # Define reference sequence
        if self.ref == '1':
            display_id = self.idA
        else:
            display_id = self.idB
        
        # Try to use logomaker
        try:
            self._display_logo_logomaker(display_df, x, x_labels, display_id)
        except ImportError as e:
            print(f"Error using logomaker for logo plot: {e}")
            print("Please install logomaker with 'pip install logomaker'")
            # Fall back to original methods if logomaker is not available
            if use_plotly.upper() in ['Y', 'YES'] and tools.get("plotly", False):
                try:
                    self._display_logo_plotly(display_df, x, x_labels, display_id)
                except ImportError as e:
                    print(f"Error using Plotly for logo plot: {e}")
                    if tools.get("matplotlib", False):
                        print("Falling back to Matplotlib...")
                        self._display_logo_matplotlib(display_df, x, x_labels, display_id)
                    else:
                        print("No visualization tools available. Please install matplotlib or plotly.")
            # Fall back to Matplotlib if available
            elif tools.get("matplotlib", False):
                try:
                    self._display_logo_matplotlib(display_df, x, x_labels, display_id)
                except ImportError as e:
                    print(f"Error using Matplotlib for logo plot: {e}")
                    if tools.get("plotly", False):
                        print("Falling back to Plotly...")
                        self._display_logo_plotly(display_df, x, x_labels, display_id)
                    else:
                        print("No visualization tools available. Please install matplotlib or plotly.")
            else:
                print("No visualization tools available. Please install matplotlib, plotly, or logomaker.")

    def _display_logo_logomaker(self, display_df: pd.DataFrame, x: pd.Series, x_labels: pd.Series, display_id: str) -> None:
        """
        Displays a sequence logo plot using the logomaker library, showing letter stacks instead of colored bars.
        
        Args:
            display_df: The filtered dataframe containing sequence data.
            x: The x-coordinates for plotting.
            x_labels: Original reference positions for x-axis labels.
            display_id: ID of reference sequence for display.
        """

        # Create figure and axes for two sequence logos
        fig, axes = plt.subplots(2, 1, figsize=(15, 8), sharex=True, 
                                gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.1})
        
        # Convert frequency dictionaries to logomaker-compatible dataframes
        logo_df_A = pd.DataFrame(index=range(len(x)))
        logo_df_B = pd.DataFrame(index=range(len(x)))
        
        # Process each position
        for pos_idx in range(len(x)):
            pos = x.iloc[pos_idx]
            pos_label = x_labels.iloc[pos_idx]
            
            # Get frequency dictionaries from dataframe
            try:
                # Check the type of the frequency data
                freq_a_data = display_df[self.FREQUENCIES_A].iloc[pos_idx]
                freq_b_data = display_df[self.FREQUENCIES_B].iloc[pos_idx]
                
                # If already a dictionary, use directly
                if isinstance(freq_a_data, dict):
                    freq_a = freq_a_data
                # If string, evaluate it
                elif isinstance(freq_a_data, str):
                    freq_a = eval(freq_a_data)
                else:
                    print(f"Unsupported frequency data type at position {pos}: {type(freq_a_data)}. Skipping.")
                    continue
                    
                # Same for freq_b
                if isinstance(freq_b_data, dict):
                    freq_b = freq_b_data
                elif isinstance(freq_b_data, str):
                    freq_b = eval(freq_b_data)
                else:
                    print(f"Unsupported frequency data type at position {pos}: {type(freq_b_data)}. Skipping.")
                    continue
                    
            except Exception as e:
                print(f"Error processing frequency data at position {pos}: {e}. Skipping.")
                continue
            
            # Add frequency data to the dataframes
            for aa, freq in freq_a.items():
                logo_df_A.loc[pos_idx, aa] = freq
                
            for aa, freq in freq_b.items():
                logo_df_B.loc[pos_idx, aa] = freq
        
        # Fill NaN values with 0
        logo_df_A = logo_df_A.fillna(0)
        logo_df_B = logo_df_B.fillna(0)
        
        # Get amino acid color scheme (can be customized)
        aa_colors = {
            'A': '#FF9966', 'C': '#FFCC00', 'D': '#CC0000', 'E': '#CC0033', 
            'F': '#00CC00', 'G': '#CC9900', 'H': '#0099CC', 'I': '#006600',
            'K': '#660099', 'L': '#009900', 'M': '#00FF00', 'N': '#CC00CC',
            'P': '#FFCC99', 'Q': '#FF00CC', 'R': '#0000CC', 'S': '#FF6600',
            'T': '#FFFF00', 'V': '#0066CC', 'W': '#660066', 'Y': '#33CC00',
            '-': '#FFFFFF', '*': '#000000', '0': '#FFFFFF'
        }

        # Create sequence logos
        logo_A = lm.Logo(logo_df_A, ax=axes[0], color_scheme=aa_colors, font_name='DejaVu Sans')
        logo_B = lm.Logo(logo_df_B, ax=axes[1], color_scheme=aa_colors, font_name='DejaVu Sans')
        
        # Configure top plot (Group A)
        axes[0].set_ylabel(f'{self.group_a}')
        axes[0].spines['bottom'].set_visible(False)
        axes[0].tick_params(labelbottom=False)
        
        # Configure bottom plot (Group B)
        axes[1].set_xlabel(f'Indices from {display_id}')
        axes[1].set_ylabel(f'{self.group_b}')
        axes[1].spines['top'].set_visible(False)
        
        # Always use the original reference positions for x-axis labels
        axes[1].set_xticks(range(len(x)))
        axes[1].set_xticklabels(x_labels, rotation=45, ha='right')
        
        if self.show.upper() in ["YES", "Y"]:
            plt.tight_layout()
            plt.show()

        # Save logoplot
        save_path = os.path.join(self.directory, "sequence_logoplot.png")
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Logoplot saved to: {save_path}")

        # For users who want to save the figure
        return fig

    def _display_logo_matplotlib(self, display_df: pd.DataFrame, x: pd.Series, x_labels: pd.Series, display_id: str) -> None:
        """
        Displays a sequence logo plot using matplotlib.
        
        Args:
            display_df: The filtered dataframe containing sequence data.
            x: The x-coordinates for plotting.
            x_labels: Original reference positions for x-axis labels.
            display_id: ID of reference sequence for display.
        """
        if plt is None:
            raise ImportError("Matplotlib is not installed. Please install matplotlib to use this feature.")
        
        # Create figure and axes
        fig, axes = plt.subplots(2, 1, figsize=(15, 10), sharex=True, 
                                gridspec_kw={'height_ratios': [1, 1], 'hspace': 0})
        
        # Get amino acid color scheme (can be customized)
        aa_colors = {
            'A': '#FF9966', 'C': '#FFCC00', 'D': '#CC0000', 'E': '#CC0033', 
            'F': '#00CC00', 'G': '#CC9900', 'H': '#0099CC', 'I': '#006600',
            'K': '#660099', 'L': '#009900', 'M': '#00FF00', 'N': '#CC00CC',
            'P': '#FFCC99', 'Q': '#FF00CC', 'R': '#0000CC', 'S': '#FF6600',
            'T': '#FFFF00', 'V': '#0066CC', 'W': '#660066', 'Y': '#33CC00',
            '-': '#CCCCCC', '*': '#000000'
        }
        
        # Create index mapping for the x-axis
        plot_indices = list(range(len(x)))
        
        # Process each position
        for pos_idx in range(len(x)):
            pos = plot_indices[pos_idx]  # Use consecutive indices for plotting
            
            # Get frequency dictionaries from dataframe
            try:
                # Check the type of the frequency data
                freq_a_data = display_df[self.FREQUENCIES_A].iloc[pos_idx]
                freq_b_data = display_df[self.FREQUENCIES_B].iloc[pos_idx]
                
                # If already a dictionary, use directly
                if isinstance(freq_a_data, dict):
                    freq_a = freq_a_data
                # If string, evaluate it
                elif isinstance(freq_a_data, str):
                    freq_a = eval(freq_a_data)
                else:
                    print(f"Unsupported frequency data type at position {pos}: {type(freq_a_data)}. Skipping.")
                    continue
                    
                # Same for freq_b
                if isinstance(freq_b_data, dict):
                    freq_b = freq_b_data
                elif isinstance(freq_b_data, str):
                    freq_b = eval(freq_b_data)
                else:
                    print(f"Unsupported frequency data type at position {pos}: {type(freq_b_data)}. Skipping.")
                    continue
                    
            except Exception as e:
                print(f"Error processing frequency data at position {pos}: {e}. Skipping.")
                continue
            
            # Sort amino acids by frequency for better visualization
            sorted_aa_a = sorted(freq_a.items(), key=lambda x: x[1], reverse=True)
            sorted_aa_b = sorted(freq_b.items(), key=lambda x: x[1], reverse=True)
            
            # Draw stacked bars for Group A (top plot)
            bottom_a = 0
            for aa, freq in sorted_aa_a:
                if freq > 0:  # Only plot if frequency > 0
                    color = aa_colors.get(aa, '#CCCCCC')  # Default gray for unknown AAs
                    bar = axes[0].bar(pos, freq, bottom=bottom_a, width=0.8, 
                                    color=color, edgecolor='black', linewidth=0.5)
                    # Add text label if bar is large enough
                    if freq >= 0.1:
                        axes[0].text(pos, bottom_a + freq/2, aa, ha='center', va='center', 
                                    color='black', fontweight='bold')
                    bottom_a += freq
            
            # Draw stacked bars for Group B (bottom plot, negative values for visual effect)
            bottom_b = 0
            for aa, freq in sorted_aa_b:
                if freq > 0:  # Only plot if frequency > 0
                    color = aa_colors.get(aa, '#CCCCCC')  # Default gray for unknown AAs
                    bar = axes[1].bar(pos, -freq, bottom=bottom_b, width=0.8, 
                                    color=color, edgecolor='black', linewidth=0.5)
                    # Add text label if bar is large enough
                    if freq >= 0.1:
                        axes[1].text(pos, bottom_b - freq/2, aa, ha='center', va='center', 
                                    color='black', fontweight='bold')
                    bottom_b -= freq
        
        # Configure top plot (Group A)
        axes[0].set_title(f'Sequence Logo Plot - {self.group_a} (top) vs {self.group_b} (bottom)')
        axes[0].set_ylabel(f'{self.group_a} Frequency')
        axes[0].set_ylim(0, 1.05)
        axes[0].spines['bottom'].set_visible(False)
        axes[0].tick_params(labelbottom=False)
        
        # Configure bottom plot (Group B)
        axes[1].set_xlabel(f'Indices from {display_id}')
        axes[1].set_ylabel(f'{self.group_b} Frequency')
        axes[1].set_ylim(-1.05, 0)
        axes[1].spines['top'].set_visible(False)
        
        # Always use the original reference positions for x-axis labels
        axes[1].set_xticks(plot_indices)
        axes[1].set_xticklabels(x_labels, rotation=45, ha='right')
        
        # Add legend for amino acid colors
        import matplotlib.patches as mpatches
        legend_patches = [mpatches.Patch(color=color, label=aa) 
                        for aa, color in aa_colors.items() 
                        if aa in set().union(*[freq_a.keys(), freq_b.keys()])]
        
        fig.legend(handles=legend_patches, loc='upper right', bbox_to_anchor=(0.99, 0.99),
                ncol=4, fontsize='small')
        
        plt.tight_layout()
        plt.subplots_adjust(right=0.85)  # Make room for the legend
        
        if self.show.upper() in ["YES" or "Y"]:
            plt.show()

        # Save Matplotlib logoplot
        save_path = os.path.join(self.directory, "matplotlib_logoplot.png")
        plt.savefig(save_path, dpi=300, bbox_inches='tight')    

    def _display_logo_plotly(self, display_df: pd.DataFrame, x: pd.Series, x_labels: pd.Series, display_id: str) -> None:
        """
        Displays an interactive sequence logo plot using Plotly.
        
        Args:
            display_df: The filtered dataframe containing sequence data.
            x: The x-coordinates for plotting.
            x_labels: Original reference positions for x-axis labels.
            display_id: ID of reference sequence for display.
        """
        if go is None:
            raise ImportError("Plotly is not installed. Please install plotly to use this feature.")
        
        # Create figure with two subplots
        fig = go.Figure()
        
        # Define amino acid color scheme (same as matplotlib version)
        aa_colors = {
            'A': '#FF9966', 'C': '#FFCC00', 'D': '#CC0000', 'E': '#CC0033', 
            'F': '#00CC00', 'G': '#CC9900', 'H': '#0099CC', 'I': '#006600',
            'K': '#660099', 'L': '#009900', 'M': '#00FF00', 'N': '#CC00CC',
            'P': '#FFCC99', 'Q': '#FF00CC', 'R': '#0000CC', 'S': '#FF6600',
            'T': '#FFFF00', 'V': '#0066CC', 'W': '#660066', 'Y': '#33CC00',
            '-': '#CCCCCC', '*': '#000000'
        }
        
        # Create index mapping for the x-axis
        plot_indices = list(range(len(x)))
        
        # Process each position
        for pos_idx in range(len(x)):
            pos = plot_indices[pos_idx]  # Use consecutive indices for plotting
            
            # Get frequency dictionaries from dataframe
            try:
                # Check the type of the frequency data
                freq_a_data = display_df[self.FREQUENCIES_A].iloc[pos_idx]
                freq_b_data = display_df[self.FREQUENCIES_B].iloc[pos_idx]
                
                # If already a dictionary, use directly
                if isinstance(freq_a_data, dict):
                    freq_a = freq_a_data
                # If string, evaluate it
                elif isinstance(freq_a_data, str):
                    freq_a = eval(freq_a_data)
                else:
                    print(f"Unsupported frequency data type at position {pos}: {type(freq_a_data)}. Skipping.")
                    continue
                    
                # Same for freq_b
                if isinstance(freq_b_data, dict):
                    freq_b = freq_b_data
                elif isinstance(freq_b_data, str):
                    freq_b = eval(freq_b_data)
                else:
                    print(f"Unsupported frequency data type at position {pos}: {type(freq_b_data)}. Skipping.")
                    continue
                    
            except Exception as e:
                print(f"Error processing frequency data at position {pos}: {e}. Skipping.")
                continue
            
            # Sort amino acids by frequency for better visualization
            sorted_aa_a = sorted(freq_a.items(), key=lambda x: x[1], reverse=True)
            sorted_aa_b = sorted(freq_b.items(), key=lambda x: x[1], reverse=True)
            
            # Draw stacked bars for Group A (positive values)
            bottom_a = 0
            for aa, freq in sorted_aa_a:
                if freq > 0:  # Only plot if frequency > 0
                    color = aa_colors.get(aa, '#CCCCCC')  # Default gray for unknown AAs
                    fig.add_trace(go.Bar(
                        x=[pos],
                        y=[freq],
                        base=[bottom_a],
                        name=aa,
                        marker_color=color,
                        marker_line=dict(width=1, color='black'),
                        text=aa,
                        textposition='inside',
                        hoverinfo='text',
                        hovertext=f"{aa}: {freq:.3f}",
                        showlegend=False if pos_idx > 0 else True,
                        legendgroup=aa
                    ))
                    bottom_a += freq
            
            # Draw stacked bars for Group B (negative values for visual effect)
            bottom_b = 0
            for aa, freq in sorted_aa_b:
                if freq > 0:  # Only plot if frequency > 0
                    color = aa_colors.get(aa, '#CCCCCC')  # Default gray for unknown AAs
                    fig.add_trace(go.Bar(
                        x=[pos],
                        y=[-freq],  # Negative value for display below x-axis
                        base=[bottom_b],
                        name=aa,
                        marker_color=color,
                        marker_line=dict(width=1, color='black'),
                        text=aa,
                        textposition='inside',
                        hoverinfo='text',
                        hovertext=f"{aa}: {freq:.3f}",
                        showlegend=False if pos_idx > 0 else True,
                        legendgroup=aa
                    ))
                    bottom_b -= freq
        
        # Configure layout
        fig.update_layout(
            title=f'Sequence Logo Plot - {self.group_a} (top) vs {self.group_b} (bottom)',
            xaxis_title=f'Indices from {display_id}',
            yaxis_title='Frequency',
            barmode='relative',
            plot_bgcolor='white',
            height=700,
            legend_title='Amino Acids',
            xaxis=dict(
                tickmode='array',
                tickvals=plot_indices,  # Always use consecutive indices for tick positions
                ticktext=list(x_labels),  # Always use original reference positions for labels
                tickangle=45
            ),
            yaxis=dict(
                range=[-1.05, 1.05],
                zeroline=True,
                zerolinewidth=2,
                zerolinecolor='black'
            )
        )
        
        # Add horizontal line at y=0
        fig.add_shape(
            type="line",
            x0=min(plot_indices)-0.5,
            y0=0,
            x1=max(plot_indices)+0.5,
            y1=0,
            line=dict(color="black", width=2)
        )
        
        # Add text annotations for group labels
        fig.add_annotation(
            x=min(plot_indices),
            y=0.9,
            text=self.group_a,
            showarrow=False,
            font=dict(size=14, color="black"),
            xanchor="left"
        )
        
        fig.add_annotation(
            x=min(plot_indices),
            y=-0.9,
            text=self.group_b,
            showarrow=False,
            font=dict(size=14, color="black"),
            xanchor="left"
        )
        
        if self.show.upper() in ["YES", "Y"]:
            fig.show()

    def display(self, use_plotly: str = 'Y', nuc_plot: str = 'Y', pm_plot: str = 'Y', create_lp: str = 'Y') -> None:
        """
        Display results using either Plotly or Matplotlib based on availability and preference.
        
        Args:
            use_plotly: Whether to prefer Plotly over Matplotlib if both are available.
            nuc_plot: Whether to display nucleotide plots when available.
            pm_plot: Whether to display point mutation plots. If None, uses self.point_mutations.
        """
        # Define reference sequence
        if self.ref == '1':
            display_id = self.idA
        else:
            display_id = self.idB
        
        # Check which display tools are available
        tools = self.check_display_tools()
        
        # Generate display dataframe once to avoid redundant processing
        display_df = self._generate_display_df()

        if create_lp.upper() in ['Y', 'YES']:
            self._create_logo_plot(display_df, use_plotly)
        
        # Use Plotly if available and preferred
        if use_plotly.upper() in ['Y', 'YES'] and tools.get("plotly", False):
            try:
                self._display_plotly(display_df, nuc_plot=nuc_plot, pm_plot=pm_plot, id=display_id)
            except ImportError as e:
                print(f"Error using Plotly: {e}")
                if tools.get("matplotlib", False):
                    print("Falling back to Matplotlib...")
                    self._display_matplotlib(display_df, pm_plot=pm_plot, id=display_id)
                else:
                    print("No visualization tools available. Please install matplotlib or plotly.")
        # Fall back to Matplotlib if available
        elif tools.get("matplotlib", False):
            try:
                self._display_matplotlib(display_df, pm_plot=pm_plot, id=display_id)
            except ImportError as e:
                print(f"Error using Matplotlib: {e}")
                if tools.get("plotly", False):
                    print("Falling back to Plotly...")
                    self._display_plotly(display_df, pm_plot=pm_plot, id=display_id)
                print("No visualization tools available. Please install matplotlib or plotly.")
        else:
            print("No visualization tools available. Please install matplotlib or plotly.")
            
#endregion

# region Input
class Input:
    #Set up argument parser for command-line parameters
    parser = argparse.ArgumentParser(description="Comparison of sequence groups in FASTA files.")

    # Required argument: first FASTA file path
    parser.add_argument("file1", type=str, help="Path to the first FASTA file")

    # Optional arguments for additional input files and parameters
    parser.add_argument("--file2", type=str, default=None, help="Path to the second file (fasta or txt)")
    parser.add_argument("--startstop", type=str, default='n', help="Translate sequence from defined start to defined stop (y/n)")
    parser.add_argument("--start", type=int, default=0, help="optional translation start")
    parser.add_argument("--stop", type=int, default=None, help="optional translation stop")
    parser.add_argument("--conservation_score", type=str, default="y", help="Calculate score with conservation or not")
    parser.add_argument("--id1", type=str, default=None, help="ID for reference sequence 1")
    parser.add_argument("--id2", type=str, default=None, help="ID for reference sequence 2")
    parser.add_argument("--top", type=int, default=20, help="Number of top scores to display")
    parser.add_argument("--lb", type=int, default=None, help="Lower bound")
    parser.add_argument("--ub", type=int, default=None, help="Upper bound")
    parser.add_argument("--positions", type=int, nargs='+', default=[], help="Specify individual sequence positions, separated by space key")
    parser.add_argument("--refseq", type=str, default="1", help="Reference sequence (1 or 2)")
    parser.add_argument("--save", type=str, default="y", help="Save results? (y/n)")
    parser.add_argument("--directory", type=str, default=os.path.expanduser("~"), help="Directory to save results")
    parser.add_argument("--show", type=str, default="N", help="Display results? (y/n)")
    parser.add_argument("--use_plotly", type=str, default="y", help="Display results with Plotly? (y/n)")
    parser.add_argument("--nuc_plot", type=str, default="n", help="Show seperate graph with triplets (Plotly only)? (y/n)")
    parser.add_argument("--pm_plot", type=str, default="n", help="Show seperate graph with point mutations? (y/n)")
    parser.add_argument("--logo_plot", type=str, default="y", help="Show amino acid logo plot? (y/n)")
    parser.add_argument("--name1", type=str, default="Group 1", help="Name of the first group")
    parser.add_argument("--name2", type=str, default="Group 2", help="Name of the second group")

# endregion

# region Main
def main():
    # Parse the command-line arguments
    args = Input.parser.parse_args()

    processor = FileProcessor()
    if args.file2:
        seq_A, seq_B = processor.process_files(args.file1, args.file2)
    else:
        seq_A, seq_B = processor.process_files(args.file1)

    type = processor.fasta_file1.type

    if type == 'NUCLEOTIDE':
        translate_A = TranslateSequence(seq_A, StartStop=args.startstop, start=args.start, stop=args.stop)
        translate_B = TranslateSequence(seq_B, StartStop=args.startstop, start=args.start, stop=args.stop)
        nuc_A, pep_A = translate_A.nuc_seq, translate_A.pep_seq
        nuc_B, pep_B = translate_B.nuc_seq, translate_B.pep_seq
        align = Alignment(pep_A, pep_B, nuc_A, nuc_B)
        pep_A, pep_B = align.pep_aligned_A, align.pep_aligned_B
        nuc_A, nuc_B = align.nuc_aligned_A, align.nuc_aligned_B
        score = Scoring(pep_A, pep_B, nuc_A, nuc_B, idA=args.id1, idB=args.id2, con_score=args.conservation_score)
        results = score.results
    if type == 'PEPTIDE':
        align = Alignment(seq_A, seq_B)
        pep_A, pep_B = align.pep_aligned_A, align.pep_aligned_B
        score = Scoring(pep_A, pep_B, idA=args.id1, idB=args.id2, con_score=args.conservation_score, ref=args.refseq)
        results = score.results
    
    output = Output(results, show=args.show.upper(), top=args.top, lb=args.lb, ub=args.ub, seqpos=tuple(args.positions), group_a=args.name1, group_b=args.name2, ref=args.refseq, idA=score.idA, idB=score.idB, directory=args.directory)
    if args.save.upper() in ["Y", "YES"]: 
        output.save()
        output.display(args.use_plotly, args.nuc_plot, args.pm_plot, args.logo_plot)

    print("Done")
        
if __name__ == "__main__":
    main()
# endregion