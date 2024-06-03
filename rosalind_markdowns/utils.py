def transcription(DNA_seq: str) -> str:
    """
    This function transcribes a DNA sequence into an RNA sequence
    -----
    Parameters:
    - DNA_seq (str): a DNA sequence string to be transcribed
    
    -----
    Returns:
    - RNA_seq (str): the transcribed RNA sequence
    """
    return DNA_seq.replace("T", "U")

def DNA_complementary(DNA_seq: str) -> str:
    """
        This function returns the complementary strand of a given DNA sequence
        -----
        Parameters:
        - DNA_seq (str): the DNA sequence string
        -----
        Returns:
        - DNA_com (str): the complementary strand
    """
    complementary = {"A" : "T", "T" : "A", "C" : "G", "G" : "C"}
    DNA_com = ""
    for nucleotide in DNA_seq:
        DNA_com += complementary[nucleotide]
    return DNA_com[::-1]

def point_mutations(s_seq: str, t_seq: str) -> int:
    """
        This function counts the number of point mutations on a sequence using
        its reference sequence
        -----
        Parameters:
        - s_seq (str): the reference DNA sequence
        - t_seq (str): the sequence to compute mutations
        -----
        Returns:
        - (int): the number of point mutations in t_seq
    """
    if (len(s_seq) != len(t_seq)):
        return 0
    else:
        point_mut = 0
        for i in range(len(s_seq)):
            if s_seq[i] != t_seq[i]:
                point_mut += 1
    return point_mut

def translation(RNA_seq: str) -> str:
    """
        This function translates the RNA sequence into a protein using the
        codon table (in data/ref folder).
        -----
        Parameters:
        - RNA_seq (str): an RNA sequence to be translated
        -----
        Returns:
        - protein (str): the encoded protein sequence
    """
    with open("../data/ref/codon.txt", "r") as file:
        text = [line.replace("\n", "").split(" ") for line in file.readlines()]
        codon_table = {}
        protein = ""
        for line in text:   
            line = [x for x in line if x != '']
            for i in range(len(line) // 2):
                codon_table[line[i*2]] = line[i*2+1]
        for i in range(len(RNA_seq) // 3):
            codon = RNA_seq[i*3:i*3+3]
            protein += codon_table[codon] if codon_table[codon] != "Stop" else ""
    return protein

from typing import List

def k_mers_enumerate(seq_str: List, k: int):
    cur = ""
    cur = recurse("", seq_str, k)
    return cur
def recurse(cur, lst, k):
    answer = ""
    if k == 0:
        return cur + "\n"
    for j in range(len(lst)):
        nextChar = lst[j]
        pr = recurse(cur + nextChar, lst, k - 1)
        answer += pr
    return answer

def get_mass(file):
    mass_dict = {}
    with open(file, "r") as f:
        text = f.read().split("\n")
        for line in text:
            line = line.split(" ")
            mass_dict[line[0]] = line[-1]
    return mass_dict

