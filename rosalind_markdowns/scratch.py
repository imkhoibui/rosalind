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
with open("data/rosalind/rosalind_kmp.txt", "r") as file:
    text = file.read().split("\n")
    text = "".join(text[1:])
    speed_up_motif(text)

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

# def shared_spliced_motif(seq_1, seq_2):
#     seq_len_1 = len(seq_1)
#     seq_len_2 = len(seq_2)
#     mat = [[""]*(seq_len_1 + 1)]*(seq_len_2 + 1)
    
#     for i in range(1, len(mat)):
#         for j in range(1, len(mat[i])):
#             print(mat[i][j-1], mat[i-1][j])
#             if seq_1[i-1] == seq_2[j-1]:
#                 mat[i][j] = mat[i-1][j-1] + seq_1[i-1]
#             else:
                
#                 mat[i][j] = common(mat[i][j-1], mat [i-1][j])
#     return mat

# def common(lst1, lst2):
#     if len(lst1) > len(lst2):
#         return lst1

# print(shared_spliced_motif("AXTGCAAAAATGTCCT", "CATGCTTTTTAGTGGAGC"))