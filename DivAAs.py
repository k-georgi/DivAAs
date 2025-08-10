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
import string
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
        ('A', 'A'): 4, ('A', 'R'): -1, ('A', 'N'): -2, ('A', 'D'): -2, ('A', 'C'): 0, ('A', 'Q'): -1, ('A', 'E'): -1, ('A', 'G'): 0, ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'L'): -1, ('A', 'K'): -1, ('A', 'M'): -1, ('A', 'F'): -2, ('A', 'P'): -1, ('A', 'S'): 1, ('A', 'T'): 0, ('A', 'W'): -3, ('A', 'Y'): -2, ('A', 'V'): 0,
        ('R', 'R'): 5, ('R', 'N'): 0, ('R', 'D'): -2, ('R', 'C'): -3, ('R', 'Q'): 1, ('R', 'E'): 0, ('R', 'G'): -2, ('R', 'H'): 0, ('R', 'I'): -3, ('R', 'L'): -2, ('R', 'K'): 2, ('R', 'M'): -1, ('R', 'F'): -3, ('R', 'P'): -2, ('R', 'S'): -1, ('R', 'T'): -1, ('R', 'W'): -3, ('R', 'Y'): -2, ('R', 'V'): -3,
        ('N', 'N'): 6, ('N', 'D'): 1, ('N', 'C'): -3, ('N', 'Q'): 0, ('N', 'E'): 0, ('N', 'G'): 0, ('N', 'H'): 1, ('N', 'I'): -3, ('N', 'L'): -3, ('N', 'K'): 0, ('N', 'M'): -2, ('N', 'F'): -3, ('N', 'P'): -2, ('N', 'S'): 1, ('N', 'T'): 0, ('N', 'W'): -4, ('N', 'Y'): -2, ('N', 'V'): -3,
        ('D', 'D'): 6, ('D', 'C'): -3, ('D', 'Q'): 0, ('D', 'E'): 2, ('D', 'G'): -1, ('D', 'H'): -1, ('D', 'I'): -3, ('D', 'L'): -4, ('D', 'K'): -1, ('D', 'M'): -3, ('D', 'F'): -3, ('D', 'P'): -1, ('D', 'S'): 0, ('D', 'T'): -1, ('D', 'W'): -4, ('D', 'Y'): -3, ('D', 'V'): -3,
        ('C', 'C'): 9, ('C', 'Q'): -3, ('C', 'E'): -4, ('C', 'G'): -3, ('C', 'H'): -3, ('C', 'I'): -1, ('C', 'L'): -1, ('C', 'K'): -3, ('C', 'M'): -1, ('C', 'F'): -2, ('C', 'P'): -3, ('C', 'S'): -1, ('C', 'T'): -1, ('C', 'W'): -2, ('C', 'Y'): -2, ('C', 'V'): -1,
        ('Q', 'Q'): 5, ('Q', 'E'): 2, ('Q', 'G'): -2, ('Q', 'H'): 0, ('Q', 'I'): -3, ('Q', 'L'): -2, ('Q', 'K'): 1, ('Q', 'M'): 0, ('Q', 'F'): -3, ('Q', 'P'): -1, ('Q', 'S'): 0, ('Q', 'T'): -1, ('Q', 'W'): -2, ('Q', 'Y'): -1, ('Q', 'V'): -2,
        ('E', 'E'): 5, ('E', 'G'): -2, ('E', 'H'): 0, ('E', 'I'): -3, ('E', 'L'): -3, ('E', 'K'): 1, ('E', 'M'): -2, ('E', 'F'): -3, ('E', 'P'): -1, ('E', 'S'): 0, ('E', 'T'): -1, ('E', 'W'): -3, ('E', 'Y'): -2, ('E', 'V'): -2,
        ('G', 'G'): 6, ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'L'): -4, ('G', 'K'): -2, ('G', 'M'): -3, ('G', 'F'): -3, ('G', 'P'): -2, ('G', 'S'): 0, ('G', 'T'): -2, ('G', 'W'): -2, ('G', 'Y'): -3, ('G', 'V'): -3,
        ('H', 'H'): 8, ('H', 'I'): -3, ('H', 'L'): -3, ('H', 'K'): -1, ('H', 'M'): -2, ('H', 'F'): -1, ('H', 'P'): -2, ('H', 'S'): -1, ('H', 'T'): -2, ('H', 'W'): -2, ('H', 'Y'): 2, ('H', 'V'): -3,
        ('I', 'I'): 4, ('I', 'L'): 2, ('I', 'K'): -3, ('I', 'M'): 1, ('I', 'F'): 0, ('I', 'P'): -3, ('I', 'S'): -2, ('I', 'T'): -1, ('I', 'W'): -3, ('I', 'Y'): -1, ('I', 'V'): 3,
        ('L', 'L'): 4, ('L', 'K'): -2, ('L', 'M'): 2, ('L', 'F'): 0, ('L', 'P'): -3, ('L', 'S'): -2, ('L', 'T'): -1, ('L', 'W'): -2, ('L', 'Y'): -1, ('L', 'V'): 1,
        ('K', 'K'): 5, ('K', 'M'): -1, ('K', 'F'): -3, ('K', 'P'): -1, ('K', 'S'): 0, ('K', 'T'): -1, ('K', 'W'): -3, ('K', 'Y'): -2, ('K', 'V'): -2,
        ('M', 'M'): 5, ('M', 'F'): 0, ('M', 'P'): -2, ('M', 'S'): -1, ('M', 'T'): -1, ('M', 'W'): -1, ('M', 'Y'): -1, ('M', 'V'): 1,
        ('F', 'F'): 6, ('F', 'P'): -4, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'W'): 1, ('F', 'Y'): 3, ('F', 'V'): -1,
        ('P', 'P'): 7, ('P', 'S'): -1, ('P', 'T'): -1, ('P', 'W'): -4, ('P', 'Y'): -3, ('P', 'V'): -2,
        ('S', 'S'): 4, ('S', 'T'): 1, ('S', 'W'): -3, ('S', 'Y'): -2, ('S', 'V'): -2,
        ('T', 'T'): 5, ('T', 'W'): -2, ('T', 'Y'): -2, ('T', 'V'): 0,
        ('W', 'W'): 11, ('W', 'Y'): 2, ('W', 'V'): -3,
        ('Y', 'Y'): 7, ('Y', 'V'): -1,
        ('V', 'V'): 4
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
    Assigns sequences to multiple groups based on different criteria.
    """
    def __init__(self):
        """Initialize empty dictionary for groups."""
        self.groups: Dict[str, Dict[str, str]] = {}
    
    def assign_by_headers(self, sequences: Dict[str, str]) -> Optional[Dict[str, Dict[str, str]]]:
        """
        Assigns sequences to groups based on the last character of their headers.
        
        Args:
            sequences: Dictionary mapping sequence headers to sequence data
            
        Returns:
            Dictionary mapping group names to dictionaries of sequences,
            or None if sequences dictionary is empty
        """
        if not sequences:
            return None
        
        # Reset groups to ensure clean assignment
        self.groups.clear()
        
        # Collect all unique group identifiers from header endings
        group_identifiers = set()
        for header in sequences.keys():
            clean_header = header.strip()
            group_identifiers.add(clean_header[-1])
        
        # Sort group identifiers for consistent naming
        sorted_identifiers = sorted(group_identifiers)
        
        # Create group name mapping
        group_mapping = {}
        for i, identifier in enumerate(sorted_identifiers):
            group_name = chr(ord('A') + i)  # A, B, C, D, ...
            group_mapping[identifier] = group_name
            self.groups[group_name] = {}
        
        # Assign sequences to groups
        for header, seq in sequences.items():
            clean_header = header.strip()
            base_header = clean_header.split()[0]
            group_identifier = clean_header[-1]
            group_name = group_mapping[group_identifier]
            
            self.groups[group_name][f"{base_header} {group_name}"] = seq
        
        return self.groups
    
    def assign_by_annotation(self, 
                            sequences: Dict[str, str], 
                            annotations: Dict[str, str]) -> Optional[Dict[str, Dict[str, str]]]:
        """
        Assigns sequences to groups based on annotation data.
        
        Args:
            sequences: Dictionary mapping sequence headers to sequence data
            annotations: Dictionary mapping sequence IDs to group annotations
            
        Returns:
            Dictionary mapping group names to dictionaries of sequences,
            or None if either input dictionary is empty
        """
        if not sequences or not annotations:
            return None
        
        # Reset groups to ensure clean assignment
        self.groups.clear()
        
        # Collect all unique group annotations
        unique_annotations = set(annotations.values())
        sorted_annotations = sorted(unique_annotations)
        
        # Create group name mapping
        group_mapping = {}
        for i, annotation in enumerate(sorted_annotations):
            group_name = chr(ord('A') + i)  # A, B, C, D, ...
            group_mapping[annotation] = group_name
            self.groups[group_name] = {}
        
        # Assign sequences to groups
        for header, seq in sequences.items():
            # Extract sequence ID from header (first word)
            seq_id = header.split()[0]
            
            if seq_id in annotations:
                annotation = annotations[seq_id]
                group_name = group_mapping[annotation]
                self.groups[group_name][f"{seq_id} {group_name}"] = seq
            else:
                # Default to first group if no annotation found
                default_group = chr(ord('A'))
                if default_group not in self.groups:
                    self.groups[default_group] = {}
                self.groups[default_group][f"{seq_id} {default_group}"] = seq
        
        return self.groups
    
    def assign_separate_files(self, *file_sequences: Dict[str, str]) -> Dict[str, Dict[str, str]]:
        """
        Assigns sequences from multiple separate files to different groups.
        
        Args:
            *file_sequences: Variable number of dictionaries, each containing sequences from a file
            
        Returns:
            Dictionary mapping group names to dictionaries of sequences
        """
        # Reset groups to ensure clean assignment
        self.groups.clear()
        
        # Assign each file to a separate group
        for i, sequences in enumerate(file_sequences):
            group_name = chr(ord('A') + i)  # A, B, C, D, ...
            self.groups[group_name] = {}
            
            for header, seq in sequences.items():
                base_header = header.split()[0]
                self.groups[group_name][f"{base_header} {group_name}"] = seq

        return self.groups
    
    def assign_by_custom_criteria(self, 
                                 sequences: Dict[str, str], 
                                 grouping_function) -> Optional[Dict[str, Dict[str, str]]]:
        """
        Assigns sequences to groups based on a custom grouping function.
        
        Args:
            sequences: Dictionary mapping sequence headers to sequence data
            grouping_function: Function that takes a header and returns a group identifier
            
        Returns:
            Dictionary mapping group names to dictionaries of sequences,
            or None if sequences dictionary is empty
        """
        if not sequences:
            return None
        
        # Reset groups to ensure clean assignment
        self.groups.clear()
        
        # Collect all unique group identifiers from the grouping function
        group_identifiers = set()
        for header in sequences.keys():
            group_identifiers.add(grouping_function(header))
        
        # Sort group identifiers for consistent naming
        sorted_identifiers = sorted(group_identifiers)
        
        # Create group name mapping
        group_mapping = {}
        for i, identifier in enumerate(sorted_identifiers):
            group_name = chr(ord('A') + i)  # A, B, C, D, ...
            group_mapping[identifier] = group_name
            self.groups[group_name] = {}
        
        # Assign sequences to groups
        for header, seq in sequences.items():
            base_header = header.split()[0]
            group_identifier = grouping_function(header)
            group_name = group_mapping[group_identifier]
            
            self.groups[group_name][f"{base_header} {group_name}"] = seq
        
        return self.groups

class AnnotationReader:
    """
    Reads and parses annotation files in various formats.
    Supports multiple groups for FASTA file processing.
    """
    def __init__(self, txt_file: str):
        """
        Initializes the AnnotationReader with a text file path.
        
        Args:
            txt_file: Path to the annotation file
        """
        self.filepath = txt_file
        self.annotation = self.read_txt()
        self.groups = self.get_unique_groups()
    
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
                    # Generate group names: A, B, C, ..., Z, AA, AB, etc.
                    group = self._generate_group_name(idx)
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

    def _generate_group_name(self, index: int) -> str:
        """
        Generates group names for large numbers of groups.
        0 -> A, 1 -> B, ..., 25 -> Z, 26 -> AA, 27 -> AB, etc.
        
        Args:
            index: Index of the group
            
        Returns:
            String representation of the group name
        """
        if index < 26:
            return chr(65 + index)  # A-Z
        else:
            # For more than 26 groups, use AA, AB, AC, etc.
            first = chr(65 + (index // 26) - 1)
            second = chr(65 + (index % 26))
            return first + second

class FileProcessor:
    """
    Processes FASTA files and assigns sequences to groups based on different criteria.
    """
    def __init__(self):
        """Initialize with a GroupAssigner."""
        self.assigner = GroupAssigner()
        self.fasta_files: Dict[int, Dict[str, str]] = {}
    
    def process_files(self, 
                     fasta: Tuple[str, ...], 
                     annotation: Optional[str] = None) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Processes one or more FASTA files and optionally an annotation file.
        
        Args:
            fasta: Tuple of paths to FASTA files (at least one required)
            annotation: Path to an annotation TXT file (optional)
        
        Returns:
            Two dictionaries (group_A, group_B) with the assigned sequences
            
        Raises:
            ValueError: If no FASTA files are provided, if both multiple FASTA files 
                       and annotation file are provided, or if file formats are invalid
        """
        if not fasta or len(fasta) == 0:
            raise ValueError("Error: At least one FASTA file must be specified!")
        
        # Check for invalid combination: multiple FASTA files + annotation file
        if len(fasta) > 1 and annotation:
            raise ValueError("Error: Cannot process multiple FASTA files together with an annotation file. Please provide either multiple FASTA files OR a single FASTA file with an annotation file.")
        
        # Load FASTA files
        fasta_sequences = []
        for i, fasta_file in enumerate(fasta):
            if not fasta_file.lower().endswith((".fasta", ".fa", ".fas")):
                raise ValueError(f"Error: File '{fasta_file}' is not a valid FASTA file (.fasta, .fa, .fas).")
            
            fasta_processor = FASTAProcessor(fasta_file)
            # fasta_sequences.append(fasta_processor.sequences)
            # self.fasta_files[i] = fasta_processor.sequences
            fasta_sequences.append(fasta_processor.sequences)
            self.fasta_files[i] = fasta_processor
        
        # Case 1: Single FASTA file with annotation
        if len(fasta) == 1 and annotation:
            if not annotation.lower().endswith(".txt"):
                raise ValueError("Error: Annotation file must be a TXT file (.txt).")
            
            annotations = AnnotationReader(annotation)
            return self.assigner.assign_by_annotation(fasta_sequences[0], annotations.annotation)
        
        # Case 2: Single FASTA file without annotation (group by headers)
        elif len(fasta) == 1 and not annotation:
            result = self.assigner.assign_by_headers(fasta_sequences[0])
            if result is None:
                return {}, {}
            return result
        
        # Case 3: Multiple FASTA files (no annotation allowed)
        elif len(fasta) > 1:
            return self.assigner.assign_separate_files(*fasta_sequences)
        
        # Fallback (should not be reached)
        return {}, {}
        

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
    A class for performing and managing pairwise and multiple sequence alignments.
    
    This class handles the alignment of protein sequences from multiple groups with optional 
    corresponding nucleotide sequences. It supports multiple alignment tools (MAFFT, Clustal Omega, MUSCLE)
    and provides functionality to map nucleotide sequences to protein alignments.
    
    Attributes:
        pep_groups (Dict[str, Dict[str, str]]): Dictionary of protein sequence groups
        nuc_groups (Dict[str, Dict[str, str]], optional): Dictionary of nucleotide sequence groups
        pep_aligned_groups (Dict[str, Dict[str, str]]): Dictionary of aligned protein sequences by group
        nuc_aligned_groups (Dict[str, Dict[str, str]], optional): Dictionary of nucleotide sequences mapped to protein alignments by group
        group_names (List[str]): List of group names
    """
    
    def __init__(self, pep_groups: Dict[str, Dict[str, str]], 
                 nuc_groups: Optional[Dict[str, Dict[str, str]]] = None):
        """
        Initialize an Alignment object and perform the alignment.
        
        Args:
            pep_groups: Dictionary where keys are group names and values are dictionaries of protein sequences
            nuc_groups: Dictionary where keys are group names and values are dictionaries of nucleotide sequences (optional)
        
        Example:
            pep_groups = {
                "group1": {"seq1": "MKQLEDKV", "seq2": "MKLLEDKV"},
                "group2": {"seq3": "MKQLEDKV", "seq4": "MKLLEDKV"},
                "group3": {"seq5": "MKQLEDKV"}
            }
        """
        if not pep_groups:
            raise ValueError("At least one group of protein sequences must be provided")
        
        self.pep_groups = pep_groups
        self.nuc_groups = nuc_groups
        self.group_names = list(pep_groups.keys())
        self.pep_aligned_groups: Dict[str, Dict[str, str]] = {}
        self.nuc_aligned_groups: Optional[Dict[str, Dict[str, str]]] = {}
        
        # Validate that nucleotide groups match protein groups if provided
        if nuc_groups:
            if set(nuc_groups.keys()) != set(pep_groups.keys()):
                raise ValueError("Nucleotide groups must have the same group names as protein groups")
            for group_name in pep_groups:
                if set(nuc_groups[group_name].keys()) != set(pep_groups[group_name].keys()):
                    raise ValueError(f"Sequence IDs in nucleotide group '{group_name}' must match protein group")
        
        # Perform protein alignment
        self._run_alignment()
        
        # Map nucleotide sequences to protein alignment if provided
        if nuc_groups:
            self.nuc_aligned_groups = {}
            for group_name in self.group_names:
                self.nuc_aligned_groups[group_name] = self._map_nucleotides_to_alignment(
                    self.nuc_groups[group_name], 
                    self.pep_aligned_groups[group_name]
                )

    def _merge_all_groups(self) -> Dict[str, str]:
        """
        Merge all protein sequence groups into a single dictionary for alignment.

        Returns:
            A single dictionary with all sequences from all groups.
        """
        merged_sequences = {}
        for group_name, sequences in self.pep_groups.items():
            for seq_id, sequence in sequences.items():
                # Add group identifier to avoid potential ID conflicts
                merged_key = f"{seq_id}_{group_name}"
                merged_sequences[merged_key] = sequence
        return merged_sequences
    
    @staticmethod
    def check_alignment_tools() -> Dict[str, bool]:
        """
        Check whether Clustal Omega, MUSCLE, or MAFFT is installed on the system.
        
        Returns:
            A dictionary with tool names as keys and booleans as values indicating if each tool is found.
        """
        tools = ["clustalo", "muscle", "mafft"]
        return {tool: shutil.which(tool) is not None for tool in tools}

    def _run_alignment(self) -> Dict[str, Dict[str, str]]:
        """
        Perform a multiple sequence alignment using an available alignment tool.
        """
        # Merge all sequence groups
        merged_sequences = self._merge_all_groups()

        # Check Alignment Tools
        tools = self.check_alignment_tools()
        
        if not any(tools.values()):
            raise RuntimeError("No alignment tool installed. Please install MAFFT, Clustal Omega, or MUSCLE.")

        # Create input sequence records with group information in the description
        input_records = []
        for group_name, group_seqs in self.pep_groups.items():
            for seq_id, seq in group_seqs.items():
                merged_id = f"{seq_id}_{group_name}"
                new_desc = f"Seq_{seq_id}_{group_name}"
                input_records.append(f">{merged_id} {new_desc}\n{seq}\n")
        
        # Create temporary input and output files for the alignment tool
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as input_file:
            input_filename = input_file.name
            input_file.write("".join(input_records))
        
        output_filename = tempfile.mktemp(suffix='.fasta')
        
        try:
            self._execute_alignment_tool(tools, input_filename, output_filename)
            
            # Process the output file to parse aligned sequences
            aligned_sequences = self._parse_alignment_output(output_filename)
            
            # Split the aligned sequences back into groups
            self.pep_aligned_groups = self._split_sequences_by_groups(aligned_sequences)
            
        except Exception as e:
            raise RuntimeError(f"Alignment failed: {e}")
        finally:
            # Clean up temporary files
            if os.path.exists(input_filename):
                os.remove(input_filename)
            if os.path.exists(output_filename):
                os.remove(output_filename)

        return self.pep_aligned_groups
    
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
    
    def _split_sequences_by_groups(self, sequences: Dict[str, Tuple[str, str]]) -> Dict[str, Dict[str, str]]:
        """
        Split the aligned sequences back into groups based on group identifiers.
        Extracts only the middle part of sequence IDs (space-separated).
        """
        grouped_sequences = {group: {} for group in self.group_names}
        
        for full_id, (description, aligned_seq) in sequences.items():
            try:
                # Extract group name and original sequence ID
                desc_parts = description.split('_')
                
                if len(desc_parts) >= 3:  # Expected format: "Seq_{seq_id}_{group_name}"
                    original_seq_id = '_'.join(desc_parts[1:-1])
                    group_name = desc_parts[-1]
                else:
                    # Fallback to splitting the ID itself
                    id_parts = full_id.split('_')
                    if len(id_parts) >= 2:
                        original_seq_id = '_'.join(id_parts[:-1])
                        group_name = id_parts[-1]
                    else:
                        print(f"Warning: Cannot parse group from {full_id}")
                        continue
                
                # Clean the sequence ID - keep only middle part if space-separated
                cleaned_id = original_seq_id.split()
                if len(cleaned_id) >= 3:  # Format like "A Identifier A"
                    cleaned_id = ' '.join(cleaned_id[1:-1])  # Keep middle parts
                elif len(cleaned_id) > 0:
                    cleaned_id = cleaned_id[0]  # Fallback to first part
                
                # Validate the group exists
                if group_name not in grouped_sequences:
                    print(f"Warning: Found sequence for unknown group '{group_name}'")
                    continue
                    
                grouped_sequences[group_name][cleaned_id] = aligned_seq
                
            except Exception as e:
                print(f"Error processing sequence {full_id}: {e}")
                continue
        
        return grouped_sequences

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
            if header not in nuc_sequences:
                continue
                
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
                    if nuc_pos + 3 <= len(nuc_seq):
                        triplet = nuc_seq[nuc_pos:nuc_pos+3]
                        nuc_mapping.append(triplet)
                        nuc_pos += 3
                    else:
                        # Handle case where nucleotide sequence is shorter than expected
                        remaining = nuc_seq[nuc_pos:]
                        nuc_mapping.append(remaining + 'N' * (3 - len(remaining)))
                        nuc_pos = len(nuc_seq)
            
            nuc_alignment[header] = ''.join(nuc_mapping)
        
        return nuc_alignment
    
    def get_group_names(self) -> List[str]:
        """
        Get the names of all groups.
        
        Returns:
            List of group names
        """
        return self.group_names.copy()
    
# endregion
    
# region Score
class Scoring:
    """
    A class for calculating similarity scores between multiple groups of aligned sequences.
    
    This class compares aligned protein sequences from multiple groups, with optional
    nucleotide sequences, and calculates position-specific scores based on amino acid
    differences, BLOSUM62 substitution matrix, hydrophobicity, and conservation rates.
    
    Attributes:
        BLOSUM62 (Dict): BLOSUM62 substitution matrix for scoring amino acid exchanges
        hydrophobicity (Dict): Dictionary mapping amino acids to their hydrophobicity values
        pep_groups (Dict[str, Dict[str, str]]): Dictionary of aligned protein sequence groups
        nuc_groups (Optional[Dict[str, Dict[str, str]]]): Dictionary of aligned nucleotide sequence groups
        ref_ids (Optional[Dict[str, str]]): Dictionary of reference sequence IDs for each group
        con_score (str): Flag indicating whether to include conservation in score calculation ("y"/"n")
        results (pd.DataFrame): DataFrame containing the calculated comparison results
    """
    
    def __init__(
        self, 
        pep_groups: Dict[str, Dict[str, str]], 
        nuc_groups: Optional[Dict[str, Dict[str, str]]] = None, 
        ref_ids: Optional[Dict[str, str]] = None, 
        con_score: str = 'y',
        ref: str = '1'
    ):
        """
        Initialize a Scoring object and calculate sequence comparison scores.
        
        Args:
            pep_groups: Dictionary of aligned protein sequence groups (group_name: sequences)
            nuc_groups: Dictionary of aligned nucleotide sequence groups (group_name: sequences) (optional)
            ref_ids: Dictionary of reference sequence IDs for each group (group_name: ref_id) (optional)
            con_score: Whether to include conservation in score calculation ("y"/"n")
        """
        self.BLOSUM62 = Constants.BLOSUM62
        self.hydrophobicity = Constants.hydrophobicity
        self.pep_groups = pep_groups
        self.nuc_groups = nuc_groups
        self.ref_ids = ref_ids or {}
        self.ref = ref
        self.con_score = con_score
        
        group_names = list(pep_groups.keys())
        if ref_ids is not None:
            if len(ref_ids) != len(group_names):
                raise ValueError("Number of reference IDs must match number of groups.")
            self.ref_ids = {group: ref_id for group, ref_id in zip(group_names, ref_ids)}
        else:
            self.ref_ids = {}

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
        - Sequence length consistency within and between groups
        - Matching reference sequence IDs if provided
        - Protein-nucleotide sequence consistency if provided
        
        Raises:
            ValueError: If any validation checks fail
        """
        if not self.pep_groups:
            raise ValueError("No sequence groups provided")
        
        # Check each group
        group_names = list(self.pep_groups.keys())
        first_group = group_names[0]
        first_len = len(next(iter(self.pep_groups[first_group].values())))
        
        # Check sequence consistency within each group and between groups
        for group_name, sequences in self.pep_groups.items():
            if not sequences:
                raise ValueError(f"Empty sequence dictionary for group {group_name}")
            
            # Check length consistency within group
            len_group = len(next(iter(sequences.values())))
            for seq in sequences.values():
                if len(seq) != len_group:
                    raise ValueError(f"Inconsistent sequence lengths within group {group_name}")
            
            # Check length consistency with first group
            if len_group != first_len:
                raise ValueError(f"Sequence length mismatch between group {first_group} ({first_len}) and group {group_name} ({len_group})")
            
            # Check reference ID if provided
            if group_name in self.ref_ids and self.ref_ids[group_name] not in sequences:
                raise ValueError(f"Reference sequence ID {self.ref_ids[group_name]} not found in group {group_name}")
        
        # Check nucleotide sequences if provided
        if self.nuc_groups:
            if set(self.nuc_groups.keys()) != set(self.pep_groups.keys()):
                raise ValueError("Nucleotide groups don't match protein groups")
            
            for group_name in self.pep_groups:
                pep_seqs = self.pep_groups[group_name]
                nuc_seqs = self.nuc_groups[group_name]
                
                # Check if keys match
                if set(pep_seqs.keys()) != set(nuc_seqs.keys()):
                    raise ValueError(f"Mismatch between protein and nucleotide sequence IDs in group {group_name}")
                
                # Check length consistency (nucleotide should be 3x protein length)
                for key in pep_seqs:
                    pep_len = len(pep_seqs[key].replace("-", ""))
                    nuc_len = len(nuc_seqs[key].replace("-", ""))
                    if nuc_len != pep_len * 3:
                        raise ValueError(f"Nucleotide sequence length for {key} in group {group_name} is not 3x the protein length")

    def _define_RefSeqs(self) -> Dict[str, str]:
        """
        Establish reference sequences based on provided IDs or use the first sequence available.
        
        Returns:
            Dictionary containing reference sequences for each group
        """
        ref_seqs = {}
        for group_name, sequences in self.pep_groups.items():
            if group_name in self.ref_ids and self.ref_ids[group_name] in sequences:
                ref_seqs[group_name] = sequences[self.ref_ids[group_name]]
            else:
                ref_seqs[group_name] = next(iter(sequences.values()))
                self.ref_ids[group_name] = next(iter(sequences.keys()))
        return ref_seqs

    def _update_score(self, aa1: str, aa2: str, score: float) -> float:
        """
        Update the score based on amino acid comparison.
        """
        if aa1 == "-" or aa2 == "-" or aa1 == aa2:
            return score
        
        blosum_score = self.BLOSUM62.get((aa1, aa2), self.BLOSUM62.get((aa2, aa1), -5))
        
        if aa1 in self.hydrophobicity and aa2 in self.hydrophobicity:
            chem_diff = abs(self.hydrophobicity[aa1] - self.hydrophobicity[aa2])
            score += abs(blosum_score) * (1 + chem_diff)
        else:
            score += abs(blosum_score)
            
        return score

    def _calculate_pairwise_score(
        self,
        seq_arrays: Dict[str, List[str]],
        most_commons: Dict[str, str],
        conservations: Dict[str, float],
    ) -> float:
        """
        Calculate pairwise scores between all groups.
        """
        total_score = 0.0
        group_names = list(seq_arrays.keys())
        
        # Compare all unique pairs of groups
        for i in range(len(group_names)):
            for j in range(i + 1, len(group_names)):
                group1 = group_names[i]
                group2 = group_names[j]
                
                score = self._calculate_single_pair_score(
                    seq_arrays[group1],
                    seq_arrays[group2],
                    most_commons[group1],
                    most_commons[group2],
                    conservations[group1],
                    conservations[group2]
                )
                
                total_score += score
        
        return total_score

    def _calculate_single_pair_score(
        self,
        seq_array_A: List[str],
        seq_array_B: List[str],
        most_common_A: str,
        most_common_B: str,
        conservation_A: float,
        conservation_B: float,
    ) -> float:
        """
        Calculate score for a single pair of groups.
        """
        if not seq_array_A or not seq_array_B:
            raise ValueError("Empty sequence arrays")

        score = 0.0
        counter_A = collections.Counter(seq_array_A)
        counter_B = collections.Counter(seq_array_B)

        # Match identical amino acids
        for aa in list(counter_A.keys()):
            if aa in counter_B:
                num_identical_pairs = min(counter_A[aa], counter_B[aa])
                for _ in range(num_identical_pairs):
                    score = self._update_score(aa, aa, score)
                counter_A[aa] -= num_identical_pairs
                counter_B[aa] -= num_identical_pairs

        # Sort remaining by hydrophobicity
        remaining_A = sorted(
            list(counter_A.elements()),
            key=lambda aa: Constants.hydrophobicity.get(aa, 0),
            reverse=True
        )
        remaining_B = sorted(
            list(counter_B.elements()),
            key=lambda aa: Constants.hydrophobicity.get(aa, 0),
            reverse=True
        )

        # Match remaining amino acids
        len_remaining_A = len(remaining_A)
        len_remaining_B = len(remaining_B)
        max_len_remaining = max(len_remaining_A, len_remaining_B)

        for j in range(max_len_remaining):
            aa1, aa2 = None, None

            if j < len_remaining_A and j < len_remaining_B:
                aa1 = remaining_A[j]
                aa2 = remaining_B[j]
            elif j < len_remaining_A:
                aa1 = remaining_A[j]
                aa2 = most_common_B
            elif j < len_remaining_B:
                aa1 = most_common_A
                aa2 = remaining_B[j]

            if aa1 is not None and aa2 is not None:
                score = self._update_score(aa1, aa2, score)

        # Adjust by conservation if needed
        if hasattr(self, 'con_score') and self.con_score.upper() in ["Y", "YES"]:
            if conservation_A is not None and conservation_B is not None:
                score = score * conservation_A * conservation_B

        return score

    def _update_Common_Conservation(
        self, 
        array: List[str],
        return_conservation: bool = True,
        return_frequencies: bool = True
    ) -> Tuple[str, Optional[float], Optional[Dict[str, float]]]:
        """
        Calculate most common element, conservation rate, and frequency dictionary.
        """
        counter = collections.Counter(array)
        max_count = max(counter.values()) if counter else 0
        
        most_common_elements = [key for key, count in counter.items() if count == max_count]
        most_common = random.choice(most_common_elements) if most_common_elements else ""
        
        conservation_rate = None
        if return_conservation:
            conservation_rate = float(max_count / len(array)) if len(array) > 0 else 0.0
        
        freq_dict = None
        if return_frequencies:
            freq_dict = {}
            total = len(array) if len(array) > 0 else 1
            for element, count in counter.items():
                freq_dict[element] = float(count / total)
                
        return most_common, conservation_rate, freq_dict

    def _remove_gaps(self, data_dict: Dict[str, Any]) -> Dict[str, Any]:
        """
        Remove positions where all groups have gaps in the reference sequences.
        Handles mixed data types (lists and strings) properly.
        """
        result_dict = data_dict.copy()
        group_names = list(self.pep_groups.keys())
        
        # Check if nucleotide data is present
        has_triplets = all(f'MostCommonTriplet_{group}' in result_dict for group in group_names)
        
        # Get reference sequences
        ref_seqs = result_dict['ref_seqs']
        seq_length = len(ref_seqs[group_names[0]])
        
        # Identify positions to keep (where at least one group has a non-gap)
        positions_to_keep = [
            i for i in range(seq_length)
            if any(ref_seqs[group][i] != "-" for group in group_names)
        ]
        
        # Filter list-type data
        for key in list(result_dict.keys()):
            if key == 'ref_seqs':
                continue  # Handle separately
                
            if isinstance(result_dict[key], dict):
                # Handle dictionary of lists (e.g., MostCommon, Conservation)
                for sub_key in list(result_dict[key].keys()):
                    if isinstance(result_dict[key][sub_key], list):
                        result_dict[key][sub_key] = [
                            result_dict[key][sub_key][i] for i in positions_to_keep
                        ]
            elif isinstance(result_dict[key], list):
                # Handle direct lists (e.g., Score)
                result_dict[key] = [result_dict[key][i] for i in positions_to_keep]
        
        # Handle nucleotide data if present
        if has_triplets:
            for group in group_names:
                triplet_key = f'MostCommonTriplet_{group}'
                result_dict[triplet_key] = [
                    result_dict[triplet_key][i] for i in positions_to_keep
                ]
            
            # Recalculate point mutations after filtering
            if len(group_names) >= 2:
                triplets_A = result_dict[f'MostCommonTriplet_{group_names[0]}']
                triplets_B = result_dict[f'MostCommonTriplet_{group_names[1]}']
                result_dict['PointMutations'] = self._triplets_with_point_mutation(triplets_A, triplets_B)
        
        # Filter reference sequences (strings)
        for group in group_names:
            result_dict['ref_seqs'][group] = ''.join(
                ref_seqs[group][i] for i in positions_to_keep
            )
        
        return result_dict

    def _define_reference_positions(self, ref_seq: str) -> List[Union[int, str]]:
        """
        Create a list of reference positions, handling gaps.
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
    
    def _triplets_with_point_mutation(self, first_triplets: List[str], second_triplets: List[str]) -> List[str]:
        """
        Check if triplet differences could be due to single point mutations.
        """
        point_mutations = []
        for i in range(len(first_triplets)):
            differences = sum(1 for a, b in zip(first_triplets[i], second_triplets[i]) if a != b)
            point_mutations.append('Y' if differences == 1 else 'N')
        return point_mutations

    def _calculate(self) -> pd.DataFrame:
        """
        Perform the comparison calculation and store the result.
        """
        self.results = self._compare_seq_groups()
        return self.results

    def _compare_seq_groups(self) -> pd.DataFrame:
        """
        Compare multiple groups of sequences and calculate position-specific scores.
        """
        group_names = list(self.pep_groups.keys())
        ref_seqs = self._define_RefSeqs()
        seq_length = len(next(iter(ref_seqs.values())))
        
        # Initialize data structures
        data_dict = {
            'ref_seqs': ref_seqs,
            'MostCommon': {group: [] for group in group_names},
            'Conservation': {group: [] for group in group_names},
            'Frequencies': {group: [] for group in group_names},
            'Score': []
        }
        
        if self.nuc_groups:
            data_dict['MostCommonTriplets'] = {group: [] for group in group_names}
        
        # Process each position
        nuc_pos = 0
        for pos in range(seq_length):
            # Collect amino acids for all groups at this position
            seq_arrays = {
                group: [seq[pos] for seq in self.pep_groups[group].values()]
                for group in group_names
            }
            
            # Update most common, conservation, and frequencies for each group
            for group in group_names:
                most_common, conservation, frequencies = self._update_Common_Conservation(seq_arrays[group])
                data_dict['MostCommon'][group].append(most_common)
                data_dict['Conservation'][group].append(conservation)
                data_dict['Frequencies'][group].append(frequencies)
            
            # Process nucleotide sequences if provided
            if self.nuc_groups:
                for group in group_names:
                    triplets = [seq[nuc_pos:nuc_pos+3] for seq in self.nuc_groups[group].values()]
                    most_common_triplet, _, _ = self._update_Common_Conservation(
                        triplets, return_conservation=False, return_frequencies=False
                    )
                    data_dict['MostCommonTriplets'][group].append(most_common_triplet)
                
                nuc_pos += 3
            
            # Calculate pairwise scores
            score = self._calculate_pairwise_score(
                seq_arrays,
                {group: data_dict['MostCommon'][group][-1] for group in group_names},
                {group: data_dict['Conservation'][group][-1] for group in group_names}
            )
            data_dict['Score'].append(score)
        
        # Remove gap positions
        cleaned_data = self._remove_gaps(data_dict)
        
        # Prepare DataFrame
        positions = list(range(1, len(cleaned_data['Score']) + 1))
        ref_positions = {
            group: self._define_reference_positions(cleaned_data['ref_seqs'][group])
            for group in group_names
        }
        
        # Normalize scores
        max_group_size = max(len(seqs) for seqs in self.pep_groups.values())
        normalized_scores = [score / max_group_size for score in cleaned_data['Score']]
        
        # Build DataFrame
        df_data = {
            "POSITION": positions,
            "SCORE": normalized_scores
        }
        
        # Add group-specific columns
        for group in group_names:
            df_data.update({
                f"REF_POS_{group}": ref_positions[group],
                f"COMMON_{group}": cleaned_data['MostCommon'][group],
                f"CONSERVATION_{group}": cleaned_data['Conservation'][group],
                f"FREQUENCIES_{group}": cleaned_data['Frequencies'][group]
            })
        
        # Add nucleotide data if available
        if self.nuc_groups and len(group_names) == 2:
            df_data.update({
                f"TRIPLET_{group_names[0]}": cleaned_data['MostCommonTriplets'][group_names[0]],
                f"TRIPLET_{group_names[1]}": cleaned_data['MostCommonTriplets'][group_names[1]],
                "POINT_MUTATIONS": cleaned_data.get('PointMutations', [])
            })
        
        return pd.DataFrame(df_data)

# endregion

# region Output
class Output:
    """
    A class for processing, saving, and visualizing comparison results between multiple sequence groups.
    
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
        groups: List[str] = ["A", "B"],  # Now accepts multiple groups
        ref: str = "1",
        ref_ids: Optional[Dict[str, str]] = None,
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
            groups: List of labels for the sequence groups.
            ref: Which reference sequence to use (index starting from "1").
            ref_ids: List of reference IDs for each group.
            directory: Output directory for saving files.
        """
        self.result_df = result_df
        self.show = show
        self.top = top
        self.lb = lb
        self.ub = ub
        self.seqpos = seqpos
        self.groups = groups
        self.num_groups = len(groups)
        self.ref = ref
        self.ref_ids = ref_ids
        self.directory = directory

        # Column names for better readability
        self.POSITION = 'POSITION'
        self.SCORE = 'SCORE'
        
        # Dynamically create column names for each group
        self.REF_POS = [f'REF_POS_{c}' for c in string.ascii_uppercase[:self.num_groups]]
        self.COMMON = [f'COMMON_{c}' for c in string.ascii_uppercase[:self.num_groups]]
        self.CONSERVATION = [f'CONSERVATION_{c}' for c in string.ascii_uppercase[:self.num_groups]]
        self.FREQUENCIES = [f'FREQUENCIES_{c}' for c in string.ascii_uppercase[:self.num_groups]]
        self.TRIPLET = [f'TRIPLET_{c}' for c in string.ascii_uppercase[:self.num_groups]]
        
        self.POINT_MUTATIONS = 'POINT_MUTATIONS'
        
    def save(self) -> None:
        """
        Saves the result dataframe as a CSV file in the specified directory.
        Handles cases where ref_ids uses group names (strings) as keys.
        """
        if not self.directory:
            raise ValueError("No directory specified for saving output")
            
        os.makedirs(self.directory, exist_ok=True)
            
        result_df_csv = self.result_df.copy()

        # Create the header with appropriate column names
        columns = {
            self.POSITION: "Position",
            self.SCORE: "Score",
        }

        # Handle reference positions - works with both numeric and string keys
        for i, group in enumerate(self.groups):
            # Get reference ID - try both numeric index and group name
            ref_id = self.ref_ids.get(i, self.ref_ids.get(group, f"Group_{group}"))
            columns[self.REF_POS[i]] = f"Positions_in_{ref_id}"

            # Add common residues, conservation and frequencies
            columns.update({
                self.COMMON[i]: f"Most_Common_{group}",
                self.CONSERVATION[i]: f"Conservation_{group}",
                self.FREQUENCIES[i]: f"AA_frequencies_{group}"
            })

        # Add triplet columns if available
        if self.TRIPLET[0] in self.result_df.columns:
            for i, group in enumerate(self.groups):
                columns[self.TRIPLET[i]] = f"Most_Common_Triplet_{group}"
            columns[self.POINT_MUTATIONS] = "Point_Mutations"

        # Rename columns for saving
        result_df_csv = result_df_csv.rename(columns=columns)
        
        # Save to CSV file
        output_path = os.path.join(self.directory, "CoSeD_output.csv")
        result_df_csv.to_csv(output_path, index=False)

        print(f"CSV saved to {output_path}")

    @staticmethod
    def check_display_tools() -> Dict[str, bool]:
        """
        Checks whether required visualization packages are installed.
        """
        tools = ["matplotlib", "plotly", "logomaker"]
        installed_tools = {tool: importlib.util.find_spec(tool) is not None for tool in tools}
        return installed_tools
    
    def _round_numerical_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Helper method to round numerical columns to two decimal places.
        """
        result = df.copy()
        result[self.SCORE] = result[self.SCORE].round(2)
        for i in range(self.num_groups):
            result[self.CONSERVATION[i]] = result[self.CONSERVATION[i]].round(2)
        return result

    def _generate_display_df(self) -> pd.DataFrame:
        """
        Generates a filtered and processed dataframe for visualization.
        """
        result_df = self.result_df

        # Apply filtering based on specific sequence positions if provided
        if self.seqpos and len(self.seqpos) > 0:
            ref_col = self.REF_POS[int(self.ref)-1]
            result_df[ref_col] = pd.to_numeric(result_df[ref_col], errors='coerce')
            display_df = result_df[result_df[ref_col].isin(self.seqpos)].copy()
            
            # Process reference position columns
            for i in range(self.num_groups):
                display_df[self.REF_POS[i]] = pd.to_numeric(display_df[self.REF_POS[i]], errors='coerce')
            
            display_df = display_df.dropna(subset=[self.REF_POS[0]])
            display_df[self.REF_POS[0]] = display_df[self.REF_POS[0]].astype(int).astype(str)
        
        # Apply position-based filtering if bounds are provided
        elif self.lb is not None and self.ub is not None:
            lb_val = max(1.0, float(self.lb))
            pos_max = result_df[self.POSITION].max()
            ub_val = min(float(self.ub), pos_max)
            
            ref_col = self.REF_POS[int(self.ref)-1]
            result_df[ref_col] = pd.to_numeric(result_df[ref_col], errors='coerce')
            display_df = result_df[(result_df[ref_col] >= lb_val) & (result_df[ref_col] <= ub_val)].copy()
            
            for i in range(self.num_groups):
                display_df[self.REF_POS[i]] = pd.to_numeric(display_df[self.REF_POS[i]], errors='coerce')
            
            display_df = display_df.dropna(subset=[self.REF_POS[0]])
            display_df[self.REF_POS[0]] = display_df[self.REF_POS[0]].astype(int).astype(str)

        else:
            top_scores_df = result_df.nlargest(self.top, self.SCORE).copy()
            ref_col = self.REF_POS[int(self.ref)-1]
            top_scores_df[ref_col] = pd.to_numeric(top_scores_df[ref_col], errors='coerce')
            display_df = top_scores_df.sort_values(by=ref_col).copy()
        
        display_df = self._round_numerical_columns(display_df)
        return display_df

    def _get_plot_data(self, display_df: pd.DataFrame) -> Tuple:
        """
        Extract data from the display dataframe needed for plotting.
        """
        ref_col = self.REF_POS[int(self.ref)-1]
        x_labels = display_df[ref_col]
        
        if self.lb is None or self.ub is None:
            x = pd.Series(range(len(display_df)))
        else:
            x = display_df[ref_col]
        
        # Extract data for all groups
        y_values = [display_df[self.CONSERVATION[i]].abs() * (1 if i % 2 == 0 else -1) for i in range(self.num_groups)]
        labels = [display_df[self.COMMON[i]] for i in range(self.num_groups)]
        score = display_df[self.SCORE]
        
        has_triplets = self.TRIPLET[0] in display_df.columns
        triplets = [display_df[self.TRIPLET[i]] if has_triplets else pd.Series([""] * len(labels[0])) for i in range(self.num_groups)]
        
        return x, y_values, labels, score, has_triplets, triplets, x_labels

    def _display_matplotlib(self, display_df: Optional[pd.DataFrame] = None, pm_plot: str = 'Y') -> None:
        """
        Displays bar plots of the result dataframe using matplotlib.
        """
        if plt is None:
            raise ImportError("Matplotlib is not installed.")
            
        if display_df is None:
            display_df = self._generate_display_df()
            
        x, y_values, labels, score, has_triplets, triplets, x_labels = self._get_plot_data(display_df)
        
        # Create subplots for each pair of groups
        num_plots = (self.num_groups + 1) // 2
        fig, axes = plt.subplots(num_plots, 1, figsize=(13, 6 * num_plots))
        if num_plots == 1:
            axes = [axes]
        
        colors = ['darkorange', 'darkgreen', 'darkblue', 'darkred', 'darkviolet', 'darkcyan']
        
        for plot_idx in range(num_plots):
            group1_idx = plot_idx * 2
            group2_idx = plot_idx * 2 + 1 if plot_idx * 2 + 1 < self.num_groups else None
            
            ax = axes[plot_idx]
            
            # Plot first group
            ax.bar(x, y_values[group1_idx], color=colors[group1_idx], label=self.groups[group1_idx])
            
            # Plot second group if exists
            if group2_idx is not None:
                ax.bar(x, y_values[group2_idx], color=colors[group2_idx], label=self.groups[group2_idx])
            
            # Annotate bars
            for i in range(len(x)):
                if group2_idx is not None:
                    # For paired groups
                    ax.text(x.iloc[i], y_values[group1_idx].iloc[i] - 0.05, 
                           f"{labels[group1_idx].iloc[i]}\n\n{triplets[group1_idx].iloc[i] if has_triplets else ''}", 
                           ha='center', va='top', fontsize=10, color='white')
                    ax.text(x.iloc[i], y_values[group2_idx].iloc[i] + 0.05, 
                           f"{labels[group2_idx].iloc[i]}\n\n{triplets[group2_idx].iloc[i] if has_triplets else ''}", 
                           ha='center', va='bottom', fontsize=10, color='white')
                else:
                    # For single group (odd number)
                    ax.text(x.iloc[i], y_values[group1_idx].iloc[i] - 0.05 if y_values[group1_idx].iloc[i] > 0 else y_values[group1_idx].iloc[i] + 0.05, 
                           f"{labels[group1_idx].iloc[i]}\n\n{triplets[group1_idx].iloc[i] if has_triplets else ''}", 
                           ha='center', va='top' if y_values[group1_idx].iloc[i] > 0 else 'bottom', fontsize=10, color='white')
                
                ax.text(x.iloc[i], 0.05, f"Score\n{score.iloc[i]}", ha='center', va='bottom', fontsize=8, color='black')
            
            # Configure plot layout
            ax.set_xlabel(f'Indices from {self.ref_ids[int(self.ref)-1]}')
            ax.set_ylabel('Conservation')
            ax.set_title(f'Conservation Plot: {self.groups[group1_idx]}' + 
                         (f' vs {self.groups[group2_idx]}' if group2_idx is not None else ''))
            
            if self.lb is None or self.ub is None:
                ax.set_xticks(x)
                ax.set_xticklabels(x_labels, rotation=45, ha='right')
            
            ax.axhline(0, color='black', linewidth=1)
            y_max = max(max(y.abs().max() for y in y_values), 1)
            ax.set_ylim(-y_max, y_max)
            ax.legend()
        
        plt.tight_layout()
        
        # Save figure
        if self.directory:
            save_path = os.path.join(self.directory, "matplotlib_plot.png")
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Matplotlib plot saved to: {save_path}")

        if self.show.upper() in ["YES", "Y"]:
            plt.show()

    def _display_plotly(self, display_df: Optional[pd.DataFrame] = None, nuc_plot: str = 'Y', pm_plot: str = 'Y') -> None:
        """
        Displays interactive bar plots of the result dataframe using Plotly.
        """
        if go is None:
            raise ImportError("Plotly is not installed.")
            
        if display_df is None:
            display_df = self._generate_display_df()
            
        x, y_values, labels, score, has_triplets, triplets, x_labels = self._get_plot_data(display_df)
        
        # Create amino acid conservation plots for each pair
        num_plots = (self.num_groups + 1) // 2
        colors = ['darkorange', 'darkgreen', 'darkblue', 'darkred', 'darkviolet', 'darkcyan']
        
        for plot_idx in range(num_plots):
            group1_idx = plot_idx * 2
            group2_idx = plot_idx * 2 + 1 if plot_idx * 2 + 1 < self.num_groups else None
            
            fig = go.Figure()
            
            # Add bars for first group
            fig.add_trace(go.Bar(
                x=x,
                y=y_values[group1_idx],
                name=self.groups[group1_idx],
                marker_color=colors[group1_idx],
                text=labels[group1_idx],
                hovertext=triplets[group1_idx] if has_triplets else labels[group1_idx],
                hoverinfo='text',
                textposition='inside',
            ))
            
            # Add bars for second group if exists
            if group2_idx is not None:
                fig.add_trace(go.Bar(
                    x=x,
                    y=y_values[group2_idx],
                    name=self.groups[group2_idx],
                    marker_color=colors[group2_idx],
                    text=labels[group2_idx],
                    hovertext=triplets[group2_idx] if has_triplets else labels[group2_idx],
                    hoverinfo='text',
                    textposition='inside',
                ))
            
            # Add score annotations
            fig.add_trace(go.Bar(
                x=x,
                y=[0.1] * len(x),
                marker_color='rgba(0,0,0,0)',
                marker_line=dict(color='white', width=0.0),
                text=score,
                textposition='outside',
                showlegend=False,
                hoverinfo='none'
            ))
            
            # Add zero line
            fig.add_trace(go.Scatter(
                x=[x.min(), x.max()],
                y=[0, 0],
                mode='lines',
                line=dict(color='black', width=1),
                showlegend=False,
                hoverinfo='none'
            ))
            
            # Customize layout
            y_max = max(max(y.abs().max() for y in y_values), 1)
            fig.update_layout(
                title=f'Conservation Plot: {self.groups[group1_idx]}' + 
                     (f' vs {self.groups[group2_idx]}' if group2_idx is not None else ''),
                xaxis_title=f'Indices from {self.ref_ids[string.ascii_uppercase[int(self.ref) - 1]]}',
                yaxis_title='Conservation',
                barmode='relative',
                showlegend=True,
                plot_bgcolor='white',
                height=600,
                yaxis=dict(range=[-y_max, y_max]),
                annotations=[
                    dict(
                        x=0.5,
                        y=y_max + 0.2,
                        xref='paper',
                        yref='y',
                        text='Scores calculated from BLOSUM62 and Hydrophobicity',
                        showarrow=False,
                        font=dict(size=12)
                    )
                ]
            )
            
            if self.lb is None or self.ub is None:
                fig.update_layout(
                    xaxis=dict(
                        tickmode='array',
                        tickvals=list(x),
                        ticktext=list(x_labels),
                        tickangle=45
                    )
                )
            
            # Save figure
            if self.directory:
                if importlib.util.find_spec("kaleido"):
                    save_path = os.path.join(self.directory, f"plotly_plot_{plot_idx+1}.png")
                    fig.write_image(save_path, width=800, height=600, scale=3)
                    print(f"Plotly plot {plot_idx+1} saved as PNG to: {save_path}")
                else:
                    html_path = os.path.join(self.directory, f"plotly_plot_{plot_idx+1}.html")
                    fig.write_html(html_path)
                    print(f"Plotly plot {plot_idx+1} saved as HTML to: {html_path}")

            if self.show.upper() in ["YES", "Y"]:
                fig.show()

    def _create_logo_plot(self, display_df: pd.DataFrame, use_plotly: str = 'Y') -> None:
        """
        Creates sequence logo plots for all groups.
        """
        tools = self.check_display_tools()
        
        if display_df is None:
            display_df = self._generate_display_df()
        
        x, _, _, _, _, _, x_labels = self._get_plot_data(display_df)
        display_id = self.ref_ids[string.ascii_uppercase[int(self.ref) - 1]]
        
        try:
            self._display_logo_logomaker(display_df, x, x_labels, display_id)
        except ImportError as e:
            print(f"Error using logomaker: {e}")
            if use_plotly.upper() in ['Y', 'YES'] and tools.get("plotly", False):
                try:
                    self._display_logo_plotly(display_df, x, x_labels, display_id)
                except ImportError as e:
                    print(f"Error using Plotly: {e}")
                    if tools.get("matplotlib", False):
                        self._display_logo_matplotlib(display_df, x, x_labels, display_id)
            elif tools.get("matplotlib", False):
                self._display_logo_matplotlib(display_df, x, x_labels, display_id)

    def _display_logo_logomaker(self, display_df: pd.DataFrame, x: pd.Series, x_labels: pd.Series, display_id: str) -> None:
        """
        Displays sequence logos using logomaker.
        """
        import logomaker as lm
        
        fig, axes = plt.subplots(self.num_groups, 1, figsize=(15, 3 * self.num_groups), 
                                sharex=True, gridspec_kw={'hspace': 0.1})
        if self.num_groups == 1:
            axes = [axes]
        
        # Process each group
        for group_idx in range(self.num_groups):
            ax = axes[group_idx]
            logo_df = pd.DataFrame(index=range(len(x)))
            
            # Process each position
            for pos_idx in range(len(x)):
                freq_data = display_df[self.FREQUENCIES[group_idx]].iloc[pos_idx]
                
                if isinstance(freq_data, dict):
                    freq = freq_data
                elif isinstance(freq_data, str):
                    freq = eval(freq_data)
                else:
                    continue
                
                for aa, f in freq.items():
                    logo_df.loc[pos_idx, aa] = f
            
            logo_df = logo_df.fillna(0)
            
            # Create logo plot
            logo = lm.Logo(logo_df, ax=ax, color_scheme='chemistry', font_name='Arial')
            
            # Configure plot
            ax.set_ylabel(self.groups[group_idx])
            if group_idx < self.num_groups - 1:
                ax.spines['bottom'].set_visible(False)
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel(f'Indices from {display_id}')
                ax.set_xticks(range(len(x)))
                ax.set_xticklabels(x_labels, rotation=45, ha='right')
        
        plt.tight_layout()
        
        if self.directory:
            save_path = os.path.join(self.directory, "sequence_logoplot.png")
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Logo plot saved to: {save_path}")

        if self.show.upper() in ["YES", "Y"]:
            plt.show()

    def display(self, use_plotly: str = 'Y', nuc_plot: str = 'Y', pm_plot: str = 'Y', create_lp: str = 'Y') -> None:
        """
        Display results using available visualization tools.
        """
        tools = self.check_display_tools()
        display_df = self._generate_display_df()

        if create_lp.upper() in ['Y', 'YES']:
            self._create_logo_plot(display_df, use_plotly)
        
        if use_plotly.upper() in ['Y', 'YES'] and tools.get("plotly", False):
            try:
                self._display_plotly(display_df, nuc_plot=nuc_plot, pm_plot=pm_plot)
            except ImportError as e:
                print(f"Error using Plotly: {e}")
                if tools.get("matplotlib", False):
                    self._display_matplotlib(display_df, pm_plot=pm_plot)
        elif tools.get("matplotlib", False):
            try:
                self._display_matplotlib(display_df, pm_plot=pm_plot)
            except ImportError as e:
                print(f"Error using Matplotlib: {e}")
            
#endregion

# region Input
class Input:
    #Set up argument parser for command-line parameters
    parser = argparse.ArgumentParser(description="Comparison of sequence groups in FASTA files.")

    # Required argument: first FASTA file path
    parser.add_argument("fasta", type=str, nargs='+', help="Path to at least one FASTA file")

    # Optional arguments for additional input files and parameters
    parser.add_argument("--annotation", type=str, default=None, help="Path to separate annotation file (txt)")
    parser.add_argument("--startstop", type=str, default='n', help="Translate sequence from defined start to defined stop (y/n)")
    parser.add_argument("--start", type=int, default=0, help="optional translation start")
    parser.add_argument("--stop", type=int, default=None, help="optional translation stop")
    parser.add_argument("--conservation_score", type=str, default="y", help="Calculate score with conservation or not")
    parser.add_argument("--ref_ids", type=str, nargs='+', help="IDs for reference sequence in each group")
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
    parser.add_argument("--group_names", type=str, nargs='+', help="Names for each group")

# endregion

# region Main
def main():
    args = Input.parser.parse_args()

    processor = FileProcessor()
    if args.annotation:
        seqs = processor.process_files(tuple(args.fasta), args.annotation)
    else:
        seqs = processor.process_files(tuple(args.fasta))

    first_fasta = processor.fasta_files[0]
    type = first_fasta.type

    if type == 'NUCLEOTIDE':
        translated = [
            TranslateSequence(seq, StartStop=args.startstop, start=args.start, stop=args.stop)
            for seq in seqs
        ]
        nuc_seqs = [t.nuc_seq for t in translated]
        pep_seqs = [t.pep_seq for t in translated]

        align = Alignment(pep_seqs, nuc_seqs)
        pep_aligned = align.pep_aligned_groups  
        nuc_aligned = align.nuc_aligned_groups

        score = Scoring(
            pep_aligned, 
            nuc_aligned,
            ref_ids = tuple(args.ref_ids) if args.ref_ids is not None else None,
            con_score=args.conservation_score,
            ref=args.refseq
        )
        results = score.results

    if type == 'PEPTIDE':
        align = Alignment(seqs)
        pep_aligned = align.pep_aligned_groups

        score = Scoring(
            pep_aligned,
            ref_ids = tuple(args.ref_ids) if args.ref_ids is not None else None,
            con_score=args.conservation_score,
            ref=args.refseq
        )
        results = score.results
    
    assigned_names = align.get_group_names()

    if (args.group_names is not None) and (len(assigned_names) != len(list(args.group_names))):
            raise ValueError("Error: Number of group names must fit the number of groups!")

    output = Output(
    results,
    show=args.show.upper(),
    top=args.top,
    lb=args.lb,
    ub=args.ub,
    seqpos=tuple(args.positions),
    groups=list(args.group_names) if args.group_names is not None else assigned_names,
    ref=args.refseq,
    ref_ids=score.ref_ids,
    directory=args.directory
    )

    if args.save.upper() in ["Y", "YES"]: 
        output.save()
        output.display(args.use_plotly, args.nuc_plot, args.pm_plot, args.logo_plot)

    print("Done")
        
if __name__ == "__main__":
    main()
# endregion