# region Import
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import re
import shutil
import subprocess
import tempfile
from collections import Counter
# endregion

# region Read Files
def read_fasta(fasta_file):
    """
    Liest eine FASTA-Datei ein und gibt ein Dictionary zurück,
    in dem die Header (ohne '>') die Keys und die Sequenzen die Values sind.
    """
    sequences = {}
    current_id = None
    try:
        with open(fasta_file, "r") as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    current_id = line[1:]
                    sequences[current_id] = ""
                elif current_id:
                    sequences[current_id] += line
    except Exception as e:
        print(f"Fehler beim Lesen der FASTA-Datei {fasta_file}: {e}")
    return sequences


def read_txt(txt_file):
    """
    Nimmt txt-Datei entgegen. 
    Entscheidet, entsprechend dem Aufbau des Inhalts, ob die Quelle ITOL oder FigTree ist
    Übergibt txt-Datei an entsprechende Funktion
    Gibt aus den Funktionen Dictionary mit Annotationen zurück
    """
    annotation = {}
    try:
        with open(txt_file, "r") as file:
            for line in file:
                if "DATA" in line.upper():
                    annotation = read_itol(txt_file)
                    return annotation
                elif "NEXUS" in line.upper():
                    annotation = read_figtree(txt_file)
                    return annotation
                elif "ANNOTATIONS" in line.upper():
                    annotation = read_annotation(txt_file)
                    return annotation
    except Exception as e:
        print(f"Fehler beim Lesen der Annotationsdatei {txt_file}: {e}")


def read_itol(itol_file):
    """
    Liest eine iTOL-Annotationstextdatei (z. B. color_annotation_AB.txt) ein und
    gibt ein Dictionary zurück, in dem die Sequenz-ID als Key und die zugehörige Gruppe (z. B. 'A' oder 'B') als Value gespeichert ist.
    
    Annahme: Die Datei enthält nach einigen Header-Zeilen ab der Zeile "DATA" 
    die Tab-getrennten Werte, z. B.:
    
    NP_0011902    range    #ff0000    A"""
    
    annotations = {}
    data_started = False
    try:
        with open(itol_file, "r") as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                if not data_started:
                    if "DATA" in line.upper():
                        data_started = True
                    continue
                # Nun nehmen wir an, dass die Zeile vier Felder enthält (getrennt durch Tab oder whitespace)
                if data_started:
                    parts = line.split()
                    # if len(parts) >= 4:
                    seq_id = parts[0]
                    group = parts[-1]  # Letztes Element sollte die Gruppe sein
                    annotations[seq_id] = group
    except Exception as e:
        print(f"Fehler beim Lesen der IOTL-Datei {itol_file}: {e}")
    return annotations

def read_figtree(fig_file):
    """
    Liest eine Annotationstextdatei ein (z. B. eine Newick-Datei mit eingebetteten Annotationen)
    und extrahiert Identifier und zugehörige Gruppenzuordnungen.

    Vorgehen:
      - Zuerst wird nach Annotationen mit !name gesucht. Wird "!name" gefunden,
        wird der in den Anführungszeichen angegebene Gruppenname übernommen.
      - Falls keine "!name"-Annotationen gefunden werden, wird stattdessen nach !color gesucht.
        Dabei wird der erste unterschiedliche Farbwert (z.B. ff0000) der Gruppe "A",
        der zweite der Gruppe "B" usw. zugeordnet.
      - Falls weder "!name" noch "!color" gefunden werden, wird "Keine Annotation übergeben" zurückgegeben.
    
    Rückgabe:
      Ein Dictionary, in dem die Keys die Identifier (z. B. "NP_001190266.1") und
      die Values die zugehörigen Gruppennamen sind.
    """
    with open(fig_file, "r") as f:
        data = f.read()

    # Versuch zuerst, Annotationen mit !name zu extrahieren:
    pattern_name = r"(?<!\))'([^']+)'\[.*?!name=\"(.*?)\""
    matches = re.findall(pattern_name, data)

    if matches:
        # Für jeden Treffer: Identifier aus einfachen Anführungszeichen und den zugehörigen Namen übernehmen
        annotations = {identifier: group for identifier, group in matches}
        return annotations

    # Falls keine !name-Annotationen gefunden wurden, versuche es mit !color:
    pattern_color = r"(?<!\))'([^']+)'\[.*?!color=#([0-9a-fA-F]{6})"
    matches_color = re.findall(pattern_color, data)

    if matches_color:
        annotations = {}
        unique_colors = []
        for identifier, color in matches_color:
            # Falls dieser Farbwert noch nicht registriert ist, hinzufügen
            if color not in unique_colors:
                unique_colors.append(color)
            idx = unique_colors.index(color)
            # Die erste gefundene Farbe (idx == 0) erhält Gruppe "A", zweite (idx == 1) Gruppe "B", usw.
            group = chr(65 + idx)  # 65 entspricht "A" im ASCII-Code
            annotations[identifier] = group
        return annotations

    # Falls weder !name noch !color gefunden wurden:
    return "Keine Annotation in figtree-Datei übergeben"

def read_annotation(annotation_file):
    """
    Nimmt Annotations-File mit Identifier und Gruppenzuweisung entgegen.
    Gibt ein Dictionary mit Identifier und Gruppenzuweisung zurück
    """
    annotations = {}
    try:
        with open(annotation_file, "r") as file:
            for line in file:
                line = line.strip()
                parts = line.split()
                seq_id = parts[0]
                group = parts[-1]  # Letztes Element sollte die Gruppe sein
                annotations[seq_id] = group
    except Exception as e:
        print(f"Fehler beim Lesen der Annotations-Datei {annotation_file}: {e}")
    return annotations


def process_files(file1, file2_or_filetxt=None):

    """
    Verarbeitet mindestens eine FASTA-Datei (file1). Optional kann eine zweite Datei übergeben werden:
    - Falls die zweite Datei eine FASTA-Datei ist, wird sie als Vergleichsdatei behandelt.
    - Falls die zweite Datei eine TXT-Datei ist, wird sie als Annotation-Datei interpretiert.
    
    Logik:
      - Ist file1 nicht vorhanden oder leer, wird ein Fehler ausgegeben.
      - Falls file2 vorhanden ist, werden file1 und file2 jeweils als Gruppe A und Gruppe B eingelesen.
      - Falls file2 fehlt:
           * Wenn filetxt vorhanden ist, wird sie ausgewertet und anhand der IDs (erster Teil des Headers)
             die Gruppenzuordnung aus der TXT-Datei übernommen.
           * Falls auch filetxt nicht vorhanden ist, wird versucht, anhand des Headers von file1 zu erkennen,
             ob eine Gruppenzuordnung (z. B. das letzte Zeichen 'A' oder 'B') vorliegt.
      - In jedem Fall werden zwei Dictionaries zurückgegeben: group_A und group_B.
    
    :param file1: Pflicht! Pfad zur ersten FASTA-Datei.
    :param file2_or_filetxt: (Optional) Entweder eine zweite FASTA-Datei oder eine TXT-Datei mit Annotationen.
    :return: Zwei Dictionaries (group_A, group_B) mit den zugehörigen Sequenzen.
    """

    # Prüfen, ob file1 existiert
    if not file1:
        raise ValueError("Fehler: Eine FASTA-Datei (file1) muss übergeben werden!")

    # Prüfen, ob file2_or_filetxt übergeben wurde
    if file2_or_filetxt:
        # Prüfen, ob es eine FASTA- oder TXT-Datei ist
        if file2_or_filetxt.lower().endswith((".fasta", ".fa", ".fas")):
            file2 = file2_or_filetxt
            filetxt = None
        elif file2_or_filetxt.lower().endswith(".txt"):
            file2 = None
            filetxt = file2_or_filetxt
        else:
            raise ValueError("Fehler: Die zweite Datei muss entweder eine FASTA- (.fasta, .fa, .fas) oder eine TXT-Datei (.txt) sein.")
    else:
        file2 = None
        filetxt = None

    # Fall 1: file2 ist vorhanden → file1 und file2 werden separat eingelesen
    if file2:
        group_A = read_fasta(file1)
        group_B = read_fasta(file2)
        return group_A, group_B

    # Fall 2: file2 nicht vorhanden → file1 einlesen
    sequences = read_fasta(file1)
    
    # Option 2a: Es wurde eine Annotationstextdatei übergeben
    if filetxt:
        annotations = read_txt(filetxt)
        group_A = {}
        group_B = {}
        for header, seq in sequences.items():
            # Annahme: Die erste "Wortgruppe" im Header entspricht der Sequenz-ID
            seq_id = header.split()[0]
            #seq_id = seq_id[1:]
            # Falls eine Annotation vorliegt, wird sie berücksichtigt
            if seq_id in annotations:
                if annotations[seq_id] == "A":
                    new_header = header + " A"
                    group_A[new_header] = seq
                elif annotations[seq_id] == "B":
                    new_header = header + " B"
                    group_B[new_header] = seq
                else: # Hier wäre die Alternative ihn einfach weg zu lassen
                    # Falls ein anderer Wert auftaucht, ordne default zu Gruppe A zu
                    group_A[header] = seq
            else: # Hier wäre die Alternative ihn einfach weg zu lassen
                # Keine Annotation für diesen Header gefunden → Default: Gruppe A
                group_A[header] = seq
        return group_A, group_B

    # Option 2b: Keine TXT-Datei – Prüfe, ob in file1 schon Annotationen vorhanden sind,
    # z. B. durch ein Stichwort oder anhand des letzten Zeichens im Header.
    group_A = {}
    group_B = {}

    for header, seq in sequences.items():
        # Wir nehmen an, dass ein gültiger Header mit Gruppenzuordnung mit 'A' oder 'B' endet.
        header_stripped = header.strip()
        if header_stripped and header_stripped[-1] in ("A", "B"): # Prüft, ob Header exisiert und A oder B die letzten Zeichen sind
            if header_stripped[-1] == "A":
                headerA = header + " A"
                group_A[headerA] = seq
            else:
                header = header + " A"
                group_B[header + " B"] = seq
        else: # Hier wäre die Alternative ihn einfach weg zu lassen
            # Falls keine Gruppenzuordnung erkannt wird, ordne default zu Gruppe A
            group_A[header] = seq

    return group_A, group_B
# endregion

# region Alignment
def merge_seq(sequences_A, sequences_B):
    """
    Führt zwei Dictionaries mit FASTA-Sequenzen in Gruppen A und B zusammen.

    :param sequences_A: Dictionary mit Sequenzen aus Gruppe A
    :param sequences_B: Dictionary mit Sequenzen aus Gruppe B
    :return: Dictionary mit den Gruppen A und B, jeweils mit den zugehörigen Sequenzen
    """
    merged_dict = {
        "A": sequences_A,  # Alle Sequenzen der ersten Datei unter "A"
        "B": sequences_B   # Alle Sequenzen der zweiten Datei unter "B"
    }
    return merged_dict

def check_alignment_tools():
    """
    Überprüft, ob Clustal Omega, MUSCLE oder MAFFT auf dem System installiert sind.
    
    Returns:
    dict: Dictionary mit {Tool-Name: True/False}, je nachdem, ob das Tool gefunden wurde.
    """
    tools = ["clustalo", "muscle", "mafft"]
    installed_tools = {tool: shutil.which(tool) is not None for tool in tools}
    
    return installed_tools

def run_alignment(merged_dict, tools):
    """
    Führt ein multiples Sequence-Alignment mit ClustalOmega durch und gibt ein Dictionary der ausgerichteten Sequenzen zurück.

    :param merged_dict: Dictionary mit den Gruppen A und B, erzeugt durch merge_seq
    :return: Dictionary mit Seq-ID als Key und Tupel (Beschreibung, Sequenz) als Value
    """
    # Erstelle Eingabe-SeqRecords mit Gruppensuffix in der Beschreibung
    input_records = []
    for group_key in ['A', 'B']:
        group_seqs = merged_dict.get(group_key, {})
        for seq_id, seq in group_seqs.items():  # Hier wird nur die Sequenz extrahiert
            new_desc = f"Seq_{seq_id}_{group_key}"  # Beschreibung mit Gruppenkennung
            input_records.append(f">{seq_id} {new_desc}\n{seq}\n")
    
    # Temporäre Eingabe- und Ausgabedateien erstellen
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as input_file:
        input_filename = input_file.name
        input_file.write("".join(input_records))  # Schreibe alle Sequenzen in die Datei
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as output_file:
        output_filename = output_file.name

    try:
        if "mafft" in tools:
            # cmd für MAFFT festlegen und MAFFT über subprocess ausführen
            cmd = [
            "mafft",  # MAFFT-Befehl
            "--auto",  # Automatische Einstellungen für das Alignment
            "--quiet", # Keine Ausgaben
            input_filename  # Eingabedatei
            ]

            with open(output_filename, "w") as outfile:
                subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
        
        elif "clustalo" in tools:
            # cmd für ClustalOmega festlegen und ClustalOmega über subprocess ausführen
            cmd = [
                "clustalo",  # ClustalOmega-Befehl
                "-i", input_filename,  # Eingabedatei
                "-o", output_filename,  # Ausgabedatei
                "--outfmt=fasta",  # Ausgabeformat
                "--force",  # Erzwinge die Ausführung
                "--auto",  # Automatische Einstellungen
                #"--verbose"  # Verbose-Modus
            ]

            subprocess.run(cmd, check=True)

        elif "muscle" in tools:
            # cmd für MUSCLE festlgene und MUSCLE über subprocess ausführen
            cmd = [
                "muscle",
                "-align", input_filename,
                "-output", output_filename,
                "-quiet"
            ]
        
            subprocess.run(cmd, check=True)
    
    except Exception as e:
        print(f"Kein Alignment-Tool installiert: {e}")
    
    # Verarbeite die Ausgabedatei
    aligned_sequences = {}
    with open(output_filename, 'r') as output_file:
        current_id = None
        current_desc = None
        current_seq = []
        for line in output_file:
            if line.startswith(">"):  # Header-Zeile
                if current_id is not None:
                    aligned_sequences[current_id] = (current_desc, "".join(current_seq))
                header_parts = line[1:].strip().split(maxsplit=1)
                current_id = header_parts[0]
                current_desc = header_parts[1] if len(header_parts) > 1 else ""
                current_seq = []
            else:  # Sequenz-Zeile
                current_seq.append(line.strip())
        if current_id is not None:  # Letzte Sequenz hinzufügen
            aligned_sequences[current_id] = (current_desc, "".join(current_seq))
    
    # Temporäre Dateien bereinigen
    os.remove(input_filename)
    os.remove(output_filename)
    
    return aligned_sequences
# endregion

# region Score
def split_seq(sequences):
    """
    Trennt die Eingabesequenzen in zwei Dictionarys, basierend auf dem letzten Zeichen der Beschreibung.

    Parameters:
    sequences (dict): Dictionary mit Sequenz-IDs als Keys und Tupeln (Beschreibung, Sequenz) als Values.

    Returns:
    dict: Zwei Dictionaries { "A": {id: seq}, "B": {id: seq} }
    """
    group_A = {}  # Dictionary für Gruppe A
    group_B = {}  # Dictionary für Gruppe B

    for seq_id, (description, sequence) in sequences.items():
        group_key = description.strip()[-1]  # Letztes Zeichen als Gruppen-ID nehmen

        if group_key == "A":
            group_A[seq_id] = sequence
        elif group_key == "B":
            group_B[seq_id] = sequence

    return group_A, group_B

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

hydrophobicity = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
    'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}

def compare_seq(seqA, seqB, idA=None, idB=None, con_score="y"):

    Reference_Positions_A =[]
    Reference_Positions_B =[]
    ConservationA = []
    ConservationB = []
    MostCommonA = []
    MostCommonB = []
    Score = []

    # Referenz-Sequenzen entsprechend ID oder als erste Seq der Gruppe anlegen
    if idA:
        ref_seqA = seqA[idA]
    else:
        ref_seqA = next(iter(seqA.values()))

    if idB:
        ref_seqB = seqB[idB]
    else:
        ref_seqB = next(iter(seqB.values()))

    for i in range(len(next(iter(seqA.values())))):
        #Positions.append(i + 1)

        # Arrays für AAs an Position i initialisieriern
        A = []
        B = []

        # Most Common Eintrag und Konservierungsgrad für erste Gruppe an Position i
        for seq in seqA.values():
            A.append(seq[i])
        counter = Counter(A)
        max_count = max(counter.values())
        most_common_list = [key for key, count in counter.items() if count == max_count]
        most_common = random.choice(most_common_list)
        MostCommonA.append(most_common) # ggf. bei Gleichstand AS aus Referenzseq nehmen
        if most_common == "-":
            ConservationA.append(0)
        else:
            ConservationA.append(float(max_count/(len(A))))

        # Most Common Eintrag und Konservierungsgrad für zweite Gruppe an Position i
        for seq in seqB.values():
            B.append(seq[i])
        counter = Counter(B)
        max_count = max(counter.values())
        most_common_list = [key for key, count in counter.items() if count == max_count]
        most_common = random.choice(most_common_list)
        MostCommonB.append(most_common) # ggf. bei Gleichstand AS aus Referenzseq nehmen
        if most_common == "-":
            ConservationB.append(0)
        else:
            ConservationB.append(float(max_count/(len(B))))

        # Überprüfen, welche Gruppe mehr Sequenzen hat
        if len(A) > len(B):
            for j in range(len(A)):
                
                if j < len(B):
                    aa1, aa2 = A[j], B[j]
                else:
                    aa1 = A[j]
                    aa2 = MostCommonB[i]

                if aa1 == "-" or aa2 == "-":
                    score = 0
            
                else:
                    score = 0
                    if aa1 == aa2:  # Keine Mutation -> Score = 0
                        score += 0
                        continue
                
                    # BLOSUM80-Wert bestimmen
                    blosum_score = BLOSUM62.get((aa1, aa2), BLOSUM62.get((aa2, aa1), -5))  # -5 als Standardwert für seltene Mutationen
                    
                    # Chemische Differenz entsprechend Hydrophobizität bestimmen
                    chem_diff = abs(hydrophobicity[aa1] - hydrophobicity[aa2])

                    # Gesamt-Score berechnen
                    score += abs(blosum_score) * (1 + chem_diff)

        elif len(A) < len(B):
            for j in range(len(B)):
                
                if j < len(A):
                    aa1, aa2 = A[j], B[j]
                else:
                    aa1 = MostCommonA[i]
                    aa2 = B[j]

                if aa1 == "-" or aa2 == "-":
                    score = 0
            
                else:
                    score = 0
                    if aa1 == aa2:  # Keine Mutation -> Score = 0
                        score += 0
                        continue
                
                    # BLOSUM80-Wert bestimmen
                    blosum_score = BLOSUM62.get((aa1, aa2), BLOSUM62.get((aa2, aa1), -5))  # -5 als Standardwert für seltene Mutationen
                    
                    # Chemische Differenz entsprechend Hydrophobizität bestimmen
                    chem_diff = abs(hydrophobicity[aa1] - hydrophobicity[aa2])

                    # Gesamt-Score berechnen
                    score += abs(blosum_score) * (1 + chem_diff)

        else:
            for j in range(len(A)):
                aa1, aa2 = A[j], B[j]

                if aa1 == "-" or aa2 == "-":
                    score = 0
                
                else:
                    score = 0
                    if aa1 == aa2:  # Keine Mutation -> Score = 0
                        score += 0
                        continue
                
                    # BLOSUM80-Wert bestimmen
                    blosum_score = BLOSUM62.get((aa1, aa2), BLOSUM62.get((aa2, aa1), -5))  # -5 als Standardwert für seltene Mutationen
                    
                    # Chemische Differenz entsprechend Hydrophobizität bestimmen
                    chem_diff = abs(hydrophobicity[aa1] - hydrophobicity[aa2])

                    # Score berechnen
                    score += abs(blosum_score) * (1 + chem_diff)
        
        if con_score.upper() in ["Y", "YES"]:
            score = score * ConservationA[i] * ConservationB[i] # Alternative Berechnung des Scores, die Konservierung mit aufnimmt

        # Gesamtscore berechnen
        Score.append(score)

    # Falls gap auf beiden Seiten, Eintrag löschen
    k = 0
    while k < len(ref_seqA):
        if MostCommonA[k] == "-" and MostCommonB[k] == "-":
            #del Positions[k]
            del MostCommonA[k]
            del ConservationA[k]
            del MostCommonB[k]
            del ConservationB[k]
            del Score[k]
            ref_seqA = ref_seqA[:k] + ref_seqA[k+1:]
            ref_seqB = ref_seqB[:k] + ref_seqB[k+1:]
        else:
            k += 1

    # Positionen insgesamt
    Positions = list(range(1, len(Score) + 1))

    # Score-Normierung
    score_max = max(Score)
    for i in range(len(Score)):
        Score[i] = Score[i]/score_max

    # Numerierungen entsprechend der Referenz-Sequenzen erstellen
    k = 0
    for i in range(len(ref_seqA)):
        if ref_seqA[i] == "-":
            Reference_Positions_A.append("nan")
        else:
            Reference_Positions_A.append(k + 1)
            k += 1
        
    k = 0
    for i in range(len(ref_seqB)):
        if ref_seqB[i] == "-":
            Reference_Positions_B.append("nan")
        else:
            Reference_Positions_B.append(k + 1)
            k += 1

    # Ergebnis-Array erstellen
    result_array = np.array([Positions, Reference_Positions_A, Reference_Positions_B, MostCommonA, ConservationA, MostCommonB, ConservationB, Score]).T

    return result_array
# endregion

# region Output
def save_score(result_array, directory, groupA = "A", groupB = "B"):
    result_array_csv = np.array(result_array, dtype=str)

    # Speichern als CSV
    np.savetxt(f"{directory}/cops_output.csv", 
            result_array_csv, delimiter=",", fmt="%s",
            header=(f"Position,Reference_Positions_{groupA},Reference_Positions_{groupB},"
            f"Most_Common_{groupA},Conservation_{groupA},Most_Common_{groupB},Conservation_{groupB},Score"), comments="")

    print("CSV gespeichert!")

def display_array(result_array, top = 20, lb = None, ub = None, groupA = "A", groupB = "B", ref = "1"):
    # Falls lb und ub übergeben werden, wird Display_Array als Bereich zwischen den bereichsgrenzen angelegt
    # Dabei bezieht sich der Bereich auf die mit ref übergebene Referenz-Sequenz
    # Zeilen, in denen "nan" in der Referenz-Sequenz auftaucht, werden vor der Ausgabe entfernt
    if lb and ub:

        lb = float(lb+1)
        ub = float(ub+1)

        pos_max = result_array[-1, 0].astype(float)

        if lb < 1:
            lb = 1
        if ub > pos_max:
            ub = pos_max

        if ref == "2":
            display_array = result_array[(result_array[:, 2].astype(float) >= lb) & (result_array[:, 2].astype(float) <= ub)]
        else:
            display_array = result_array[(result_array[:, 1].astype(float) >= lb) & (result_array[:, 1].astype(float) <= ub)]

         # Versuche, alle Werte in Float zu konvertieren, ersetze Fehler durch NaN
        display_array[:, 1] = np.array([float(x) if str(x).replace('.', '', 1).isdigit() else np.nan for x in display_array[:, 1]])

        # Entferne Zeilen mit NaN
        display_array = display_array[~np.isnan(display_array[:, 1].astype(float))]

        # 1. Spalte 1 in Strings umwandeln
        display_array[:, 1] = display_array[:, 1].astype(str)

        # 2. Letzte zwei Zeichen entfernen (falls der String mindestens zwei Zeichen hat)
        display_array[:, 1] = np.array([x[:-2] if len(x) > 2 else x for x in display_array[:, 1]])

        # 3. In Integer umwandeln
        display_array[:, 1] = display_array[:, 1].astype(int)

    # Ohne lb und ub werden die top-Scores (default = 20) als Display-Array angelegt
    else:
        display_array = result_array[result_array[:, 7].argsort()[::-1]]
        display_array = display_array[:top]

    # Positionen welcher Ref-Seq werden als x-Achse gewählt
    if ref == "2":
        x = display_array[:,2]
    else:
        x = display_array[:,1]

    # Runden auf zwei Nachkommastellen
    display_array[:, 4] = np.round(display_array[:, 4].astype(float), 2).astype(float) # Conservation A
    display_array[:, 6] = np.round(display_array[:, 6].astype(float), 2).astype(float) # Conservation B
    display_array[:, 7] = np.round(display_array[:, 7].astype(float), 2).astype(float) # Score

    y_up = np.abs(display_array[:, 4].astype(float)) # Konservierungsgrad A
    y_down = 0-np.abs(display_array[:, 6].astype(float)) #Konservierungsgrad B

    labelsA = display_array[:, 3]
    labelsB = display_array[:, 5]
    score = display_array[:, 7]

    # Erstellen des Plots
    fig, ax = plt.subplots(figsize=(13, 8))

    # Positive Balken (Standardfarbe Blau)
    ax.bar(x, y_up, color='blue', label=groupA)

    # Negative Balken (Standardfarbe Rot)
    ax.bar(x, y_down, color='red', label=groupB)

    # Beschriftung der Balken mit den häufigsten AS
    for i in range(len(x)):
        ax.text(x[i], y_up[i] - 0.05, str(labelsA[i]), ha='center', va='top', fontsize=12, color='black')
        ax.text(x[i], y_down[i] + 0.05, str(labelsB[i]), ha='center', va='bottom', fontsize=12, color='black')
        ax.text(x[i], 0.05, str(score[i]), ha='center', va='bottom', fontsize=8, color='black')

    # Achsenbeschriftungen
    ax.set_xlabel(f'Indizes aus Gruppe {ref}')
    ax.set_ylabel('Konservierungsgrad')
    ax.set_title('Barplot der Gruppen')

    # Null-Linie für bessere Sichtbarkeit
    ax.axhline(0, color='black', linewidth=1)

    ax.set_ylim(-1, 1)

    # Legende
    ax.legend()

    plt.subplots_adjust(bottom=0.2)  # Vergrößert den unteren Rand des Plots

    # Plot anzeigen
    plt.show()
# endregion

# region Main   
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Vergleich von Sequenzgruppen in FASTA-Dateien.")

    # Erforderliches Argument
    parser.add_argument("fasta_file1", type=str, help="Pfad zur ersten FASTA-Datei")

    # Optionale Argumente
    parser.add_argument("--fasta_file2", type=str, default=None, help="Pfad zur zweiten FASTA-Datei (optional)")
    parser.add_argument("--txt_file", type=str, default=None, help="Pfad zur Annotationstextdatei (optional)")
    parser.add_argument("--conservation_score", type=str, default="y", help="Score-Berechnung mit Konservierungsgrad oder nicht")
    parser.add_argument("--id1", type=str, default=None, help="ID für Referenz-Sequenz 1")
    parser.add_argument("--id2", type=str, default=None, help="ID für Referenz-Sequenz 2")
    parser.add_argument("--top", type=int, default=20, help="Anzahl der besten Werte ausgeben")
    parser.add_argument("--lb", type=int, default=None, help="Untere Grenze")
    parser.add_argument("--ub", type=int, default=None, help="Obere Grenze")
    parser.add_argument("--refseq", type=str, default="1", help="Referenzsequenz (1 oder 2)")
    parser.add_argument("--save", type=str, default="y", help="Ergebnisse speichern? (y/n)")
    parser.add_argument("--display", type=str, default="y", help="Ergebnisse anzeigen? (y/n)")
    parser.add_argument("--directory", type=str, default=os.path.expanduser("~"), help="Speicherort für Ergebnisse")
    parser.add_argument("--name1", type=str, default="Gruppe 1", help="Name der ersten Gruppe")
    parser.add_argument("--name2", type=str, default="Gruppe 2", help="Name der zweiten Gruppe")

    args = parser.parse_args()

    # Aufruf der Funktionen mit den übergebenen Argumenten
    if args.txt_file:
        seqA, seqB = process_files(args.fasta_file1, args.txt_file)
    elif args.fasta_file2:
        seqA, seqB = process_files(args.fasta_file1, args.fasta_file2)
    else:
        seqA, seqB = process_files(args.fasta_file1)   
    merged_seq = merge_seq(seqA, seqB)
    tools = check_alignment_tools()
    aligned_seq = run_alignment(merged_seq, tools)
    alignA, alignB = split_seq(aligned_seq)
    result_array = compare_seq(alignA, alignB, idA=args.id1, idB=args.id2, con_score=args.conservation_score)

    if args.save.upper() in ["Y", "YES"]:
        save_score(result_array, args.directory, groupA=args.name1, groupB=args.name2)

    if args.display.upper() in ["Y", "YES"]:
        display_array(result_array, lb=args.lb, ub=args.ub)

# endregion