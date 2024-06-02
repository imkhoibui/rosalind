from typing import List

def insertion_sort(ls: List) -> List:
    new_ls = [ls[0]]
    for i in range(1, len(ls)):
        pivot = ls[i]
        
        if pivot >= new_ls[-1]:
            new_ls.append(pivot)
        else:
            for j in range(len(new_ls) -1, -1, -1):
                if pivot >= new_ls[j]:
                    new_ls.insert(j+1, pivot)
                    break
                elif j == 0:
                    new_ls.insert(0, pivot)
                elif pivot < new_ls[j]:
                    continue
    return new_ls

a_list = [55, 17, 23, 1, 35, 89, 43, 4, 22, 11]
print(insertion_sort(a_list))