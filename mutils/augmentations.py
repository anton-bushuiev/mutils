import numpy as np
import pandas as pd
import itertools
import warnings
import os
import math


def reverse_permutation(df, complex="#Pdb", mutstr="mutstr", ddG="ddG", num_muts="num_muts"):
    df = df.drop_duplicates(subset=[complex, mutstr], ignore_index=True)
    
    dups = len(df[df.duplicated([complex, mutstr])])
    if dups > 0:
        warnings.warn(str(dups) + " duplicated entries")
    
    reversed_mutations = []
    complexes = []
    rev_signed = []
    nums = []

    for mut in np.array(df[mutstr]):
        mutats = mut.split(",", -1)
        reverse = []
        
        for each in mutats:
            swapped = each[-1] + each[1:-1] + each[0]
            reverse.append(swapped)

        r = ','.join(reverse)
        reversed_mutations.append(r)
        nums.append(len(reverse))
        
    complexes = list(df[complex])
    rev_signed = [-x if x != 0 else x for x in df[ddG]]
    
    lists = list(zip(complexes, reversed_mutations, nums, rev_signed))
    res = pd.DataFrame(lists, columns=[complex, mutstr, num_muts, ddG])
    res["original"] = [False for x in range(len(res))]

    return res


def permutations(df, complex="#Pdb", mutstr="mutstr", ddG="ddG", num_muts="num_muts"):
    df = df.drop_duplicates(subset=[complex, mutstr], ignore_index=True)

    dups = len(df[df.duplicated([complex, mutstr])])
    if dups > 0:
        warnings.warn(str(dups) + " duplicated entries")

    sorted = df[df[num_muts]== 1].sort_values(by=[mutstr, complex], ignore_index=True, key=lambda x: x.str[0:-1])
    result = pd.DataFrame()
    char_idx = 0
    code_idx = 1

    while(char_idx < len(sorted[mutstr])):
        same_code = []
        same_code.append([sorted[mutstr][char_idx], sorted[ddG][char_idx]])
        while(code_idx < len(sorted[mutstr]) and sorted[mutstr][char_idx][1:-1] == sorted[mutstr][code_idx][1:-1] ):
            if sorted[mutstr][char_idx][0] != sorted[mutstr][code_idx][0] or sorted[complex][char_idx] != sorted[complex][code_idx] :
                break

            same_code.append([sorted[mutstr][code_idx], sorted[ddG][code_idx]])
            code_idx += 1
        
        res = proccess_permutation(same_code, sorted[complex][char_idx], complex, mutstr, ddG, num_muts)
        result = pd.concat([result, res], ignore_index=True)

        char_idx = code_idx
        code_idx += 1
    
    result["original"] = [False for x in range(len(result))]
    return result


def proccess_permutation(mutations, cmplx, complex, mutstr, ddG, num_muts):
    if len(mutations) <= 1:
        return None

    chars = {mutations[i][0][-1]:mutations[i][1] for i in range(len(mutations))}
    code = mutations[0][0][1:-1]
    perms = list(itertools.permutations(chars, 2))

    muts = []
    complexes = []
    ddGs = []
    for each in perms:
        mstr = each[0] + code + each[1]
        deltaG = round(chars[each[1]] - chars[each[0]], 6)
        muts.append(mstr)
        complexes.append(cmplx)
        ddGs.append(deltaG)
    
    nums = [1 for x in range(len(muts))]
    lists = list(zip(complexes, muts, nums, ddGs))
    res = pd.DataFrame(lists, columns=[complex, mutstr, num_muts, ddG])
    return res


def mutation_distribution(df):
    to_index = {"A": 0, "R": 1, "N": 2, "D": 3, "C": 4, 
                "Q": 5, "E": 6, "G": 7, "H": 8, "I": 9, 
                "L": 10, "K": 11, "M": 12, "F": 13, "P": 14, 
                "S": 15, "T": 16, "W": 17, "Y": 18, "V": 19}
    distribution = np.zeros((20, 20), dtype=int)
    for mut in df.mutstr:
        fro = mut[0]
        to = mut[-1]
        distribution[to_index[fro]][to_index[to]] += 1

    return distribution


def save_csvs(filepath, name, df, complex="#Pdb", mutstr="mutstr", ddG="ddG", num_muts="num_muts"):
    if num_muts not in df.columns:
        df[num_muts] = df[mutstr].str.count(',') + 1

    reverse = reverse_permutation(df, complex, mutstr, ddG, num_muts)
    permute = permutations(df, complex, mutstr, ddG, num_muts)

    df_reverse = pd.concat([df, reverse], ignore_index=True)
    df_permute = pd.concat([df, permute], ignore_index=True)

    reverse_permute = permutations(df_reverse, complex, mutstr, ddG, num_muts)
    permute_reverse = reverse_permutation(df_permute, complex, mutstr, ddG, num_muts)

    df_reverse_permute = pd.concat([df, reverse_permute], ignore_index=True)
    df_permute_reverse = pd.concat([df, permute_reverse], ignore_index=True)
    
    os.makedirs(filepath, exist_ok=True)
    df_reverse.to_csv(filepath + "/" + name + "_reverse.csv")
    df_permute.to_csv(filepath + "/" + name + "_permute.csv")
    df_reverse_permute.to_csv(filepath + "/" + name + "_reverse_permute.csv")
    df_permute_reverse.to_csv(filepath + "/" + name + "_permute_reverse.csv")
    

def create_augmentations(df, complex="#Pdb", mutstr="mutstr", ddG="ddG", num_muts="num_muts"):
    reverse = reverse_permutation(df, complex, mutstr, ddG, num_muts)
    permute = permutations(df, complex, mutstr, ddG, num_muts)

    df_reverse = pd.concat([df, reverse], ignore_index=True)
    df_permute = pd.concat([df, permute], ignore_index=True)

    reverse_permute = permutations(df_reverse, complex, mutstr, ddG, num_muts)
    permute_reverse = reverse_permutation(df_permute, complex, mutstr, ddG, num_muts)

    df_reverse_permute = pd.concat([df, reverse_permute], ignore_index=True)
    df_permute_reverse = pd.concat([df, permute_reverse], ignore_index=True)

    dfs = [df_reverse, df_permute, df_reverse_permute, df_permute_reverse]
    lenghts = [
        ["reverse", len(df_reverse)],
        ["permute", len(df_permute)],
        ["reverse_permute", len(df_reverse_permute)],
        ["permute_reverse", len(df_permute_reverse)]
    ]
    return dfs, lenghts


def preprocess_skempi(df):# #Pdb, Mutation(s)_cleaned, Affinity_mut_parsed, Affinity_wt_parsed, Temperature
    complexes = []
    mutstrings = []
    numbers_mut = []
    ddGs = []
    originality = [True for x in range(len(df))]
    for i, row in df.iterrows():
        complexes.append(row["#Pdb"])
        mutstrings.append(row["Mutation(s)_cleaned"])
        ddGs.append(row["ddG"])
        num = len(row["Mutation(s)_cleaned"].split(",", -1))
        numbers_mut.append(num)

    lists = list(zip(complexes, mutstrings, numbers_mut, ddGs, originality))
    res = pd.DataFrame(lists, columns=["#Pdb", "Mutation(s)_cleaned", "num_muts", "ddG", "original"])
    return res
