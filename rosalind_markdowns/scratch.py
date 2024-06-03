def speed_up_motif(seq: str):
    ans = ""

    ## construct 0 array
    arr = []

    for i in range(len(seq)):
        arr.append(0)

    for i in range(1, len(seq)):
        cumm = 0
        for j in range(len(seq)):
            if j + i < len(seq) and seq[j] == seq[j+i]:
                cumm += 1
                if cumm > arr[j+i]:
                    arr[j+i] = cumm
            else:
                break
    
    with open("data/results/results_motif.txt", "w") as f:
        ans = ""
        for each in arr:
            ans += str(each) + " "
        f.write(ans)
    return ans
# with open("data/rosalind/rosalind_kmp.txt", "r") as file:
#     text = file.read().split("\n")
#     text = "".join(text[1:])
#     speed_up_motif(text)

# print(speed_up_motif("AABABCABCDABCDEABCDEF"))

def enumerate_oriented_genes(n: int):
    lst = []
    for i in range(n):
        lst.append(i + 1)
    print("Before enumerating", lst)
    return enumerate_helper(lst, "", n)

def enumerate_helper(lst, seq, k):
    ans = ""
    ## base case, stops when k == 0
    if k == 0:
        return seq + "\n"
    
    for i in range(len(lst)):
        print(lst)
        char = lst.pop(0)
        seq_pos = seq + str(char) + " "
        ans += enumerate_helper(lst, seq_pos, k - 1)
        
        neg_char = char*-1
        seq_neg = seq + str(neg_char) + " "
        ans += enumerate_helper(lst, seq_neg, k - 1)

        lst.append(char)

    return ans

# print(enumerate_oriented_genes(2))

def random_motifs(N, x, s):
    P = 1.0
    for char in s:
        if char in 'AT':
            P *= (1 - x)/2
        elif char in 'GC':
            P *= x/2
    prob = 1 - (1 - P)**N
    return prob


# print(random_motifs(82568, 0.539812, "AGCTTGCC"))


# def shared_spliced_motif(seq_1, seq_2):
#     ss_motif = ""
#     seq_len = len(seq_1)
#     score_arr = []

#     ## constructing score array
#     for i in range(seq_len):
#         score_arr.append(0)

#     ## checking for similarity
#     start = 0
#     for i in range(seq_len):
#         if seq_1[i] == seq_2[i]:
#             ss_motif += longest_common_subsequence(seq_1[start:i], seq_2[start:i])
#             ss_motif += seq_1[i]
#         else:
#             start = i
#     return

# def longest_common_subsequence(seq_1, seq_2):

#     return

def expected_restriction(n, substr, arr):
    str_len = len(substr)
    n_select = n - str_len + 1
    expected = []

    for gc_content in arr:
        GC_prob = gc_content / 2
        AT_prob = 0.5 - GC_prob
        prob = 1.0
        for char in substr:
            if char in "GC":
                prob *= GC_prob
            elif char in "AT":
                prob *= AT_prob
        expected.append(n_select*prob)

    return expected

# with open("data/rosalind/rosalind_eval.txt", "r") as file:
#     text = file.read().split("\n")
#     n = int(text[0])
#     substr = text[1]
#     arr = []
#     for char in text[2].split():
#         arr.append(float(char))

#     print(expected_restriction(n, substr, arr))

# def catalan_number(seq):

#     ## base case
#     if len(seq) == 2:
#         return 


# def alternative_splicing(n, k):
#     sum_com = 0
#     for i in range(k, n + 1):
#         comm = 1
#         for j in range(1, i + 1):
#             con = ((n - j + 1) / (j)) 
#             comm *= con
#         sum_com += comm % 1000000
#     return sum_com 

# print(alternative_splicing(169, 123))

def shared_spliced_motif(seq_1, seq_2):
    seq_len_1 = len(seq_1)
    seq_len_2 = len(seq_2)
    mat = [[""]*(seq_len_1 + 1)]*(seq_len_2 + 1)
    
    for i in range(1, len(mat)):
        for j in range(1, len(mat[i])):
            print(mat[i][j-1], mat[i-1][j])
            if seq_1[i-1] == seq_2[j-1]:
                mat[i][j] = mat[i-1][j-1] + seq_1[i-1]
            else:
                
                mat[i][j] = common(mat[i][j-1], mat [i-1][j])
    ans = traceback(mat)
    return ans

def traceback(mat):

    return mat

def common(lst1, lst2):
    if len(lst1) > len(lst2):
        return lst1

# print(shared_spliced_motif("AXTGCAAAAATGTCCT", "CATGCTTTTTAGTGGAGC"))

def set_operations(n, set_1, set_2):
    ans = ""
    full_set = set()
    for value in range(int(n)):
        full_set.add(value+1)

    ## Convert the strings to sets
    seq_1 = set_1.strip("{}")
    seq_2 = set_2.strip("{}")

    seq_1 = seq_1.split(",")
    seq_2 = seq_2.split(", ")

    set_1 = {int(element) for element in seq_1}
    set_2 = {int(element) for element in seq_2}

    ## Union
    ans += str(set_1.union(set_2)) + "\n"
    ## Intersection
    ans += str(set_1.intersection(set_2)) + "\n"

    ## Set difference 1 - 2
    ans += str(set_1 - set_2) + "\n"
    ## Set difference 2 - 1
    ans += str(set_2 - set_1) + "\n"

    ## Complementary A
    ans += str(full_set - set_1) + "\n"
    ## Complementary B
    ans += str(full_set - set_2) + "\n"
    
    with open("data/results/results_seto.txt", "w") as f:
        f.write(ans)

# with open("data/rosalind/rosalind_seto.txt", "r") as file:
#     text = file.read().split("\n")
#     n = text[0]
#     set_1 = text[1]
#     set_2 = text[2]
#     set_operations(n, set_1, set_2)

from typing import List
from utils import get_mass

def protein_infer(L_list: List):
    mass_dict = get_mass("data/ref/mass.txt")
    protein_str = ""
    for i in range(1, len(L_list)):
        residue = round(float(L_list[i]) - float(L_list[i-1]), 4)
        for key, value in mass_dict.items():
            
            if round(float(value), 4) == residue:
                protein_str += key
                break
    return protein_str

# with open("data/rosalind/rosalind_spec.txt", "r") as f:
#     text = f.read()
#     L_list = text.split("\n")[:-1]
#     print(protein_infer(L_list))

def overlap(seq_1, seq_2):
    max_overlap = 0

    min_len = min(len(seq_1), len(seq_2))
    for i in range(1, min_len + 1):
        if seq_1[-i:] == seq_2[:i]:
            if i > max_overlap:
                max_overlap = i
    
    return max_overlap
    

def genome_assembly(seqs):
    while len(seqs) > 1:
        max_overlap = 0
        best_pair = (0, 0)
        superstring = ""

        for i in range(len(seqs)):
            for j in range(len(seqs)):
                if i != j:
                    overlap_length = overlap(seqs[i], seqs[j])
                    if overlap_length > max_overlap:
                        max_overlap = overlap_length
                        best_pair = (i, j)
                        superstring = seqs[i] + seqs[j][overlap_length:] 
        i, j = best_pair
        seqs[i] = superstring
        seqs.pop(j)

    return seqs[0]

# with open("data/rosalind/rosalind_superstring.txt", "r") as file:
#     text = file.read().split(">")
#     seqs = []
#     for seq in text[1:]:
#         seq = seq.split("\n")[1:]
#         seqs.append("".join(seq))

#     print(genome_assembly(seqs))

from utils import point_mutations, DNA_complementary

def error_corrections(seqs):
    seqs_dict = {}
    correct_dict = {}
    incorrect_dict = {}

    for seq in seqs:
        complementary = DNA_complementary(seq)
        if seq in seqs_dict:
            seqs_dict[seq] += 1
        elif complementary in seqs_dict:
            seqs_dict[complementary] += 1
        elif seq not in seqs_dict:
            seqs_dict[seq] = 1

    for seq, value in seqs_dict.items():
        if value >= 2:
            correct_dict[seq] = value
        else:
            incorrect_dict[seq] = value

    ans = ""
    for seq in incorrect_dict.keys():
        for correct_seq in correct_dict.keys():
            if point_mutations(seq, correct_seq) == 1:
                ans += seq + "->" + correct_seq + "\n"
            elif point_mutations(DNA_complementary(seq), correct_seq) == 1:
                ans += seq + "->" + correct_seq + "\n"

    return ans

# with open("data/rosalind/rosalind_corr.txt", "r") as file:
#     text = file.read().split(">")[1:]
#     seqs = ["".join(seq.split("\n")[1:])for seq in text]
#     print(error_corrections(seqs))

import math
def catalan_number(seqs):
    AU_count = 0
    GC_count = 0

    for char in seqs:
        if char in "A":
            AU_count += 1
        elif char in "G":
            GC_count += 1
    
    catalan = []
    catalan.append(1)
    catalan.append(1)
    edge_num = AU_count * GC_count

    for i in range(2, edge_num):
        catalan.append(0)
        for j in range(i):
            catalan[i] += (catalan[j] * catalan[i - j - 1]) % 1000000
            if catalan[i] >= 1000000:
                catalan[i] -= 1000000
    return catalan

# with open("data/rosalind/rosalind_cat.txt", "r") as file:
#     text = file.read().split("\n")
#     text = "".join(text[1:])
#     print(catalan_number("CGGCUGCUACGCGUAAGCCGGCUGCUACGCGUAAGC"))

def comparing_spectra(seq_1, seq_2):
    difference = {}
    for num_1 in seq_1:
        for num_2 in seq_2:
            diff = float(num_1) - float(num_2)
            diff = round(diff, 5)
            if diff not in difference:
                difference[diff] = 1
            else:
                difference[diff] += 1
    multiplicity = max(difference, key=difference.get)
    max_value = difference[multiplicity]
    return max_value, multiplicity

with open("data/rosalind/rosalind_conv.txt", "r") as file:
    text = file.read().split("\n")
    seq_1 = [float(num) for num in text[0].split(" ")]
    seq_2 = [float(num) for num in text[1].split(" ")]

    print(comparing_spectra(seq_1, seq_2))
