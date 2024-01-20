import pytest
import pandas as pd
from mutils.augmentations import reverse_permutation, permutations


def test_reverse_permutation():
    data = [
        ["1A22_A_B", "DA11A,MA14V,HA18Q,RA19H,FA25A,QA29K,EA33R", 7, 0.747160, True],
        ["1MLC_AB_E", "TB58D", 1, -0.559933, True],
    ]
    expected_reversed = [
        ["1A22_A_B", "AA11D,VA14M,QA18H,HA19R,AA25F,KA29Q,RA33E", 7, -0.747160, False],
        ["1MLC_AB_E", "DB58T", 1, 0.559933, False],
    ]
    columns = ["complex", "mutstr", "num_muts", "ddG", "original"]

    df = pd.DataFrame(data, columns=columns)
    df_exp = pd.DataFrame(expected_reversed, columns=columns)

    df_rev = reverse_permutation(df, "complex", "mutstr", "ddG", "num_muts")

    assert df_rev.equals(df_exp)

def test_permutations():
    data = [
        ["1A22_A_B", "DA11A,MA14V,HA18Q,RA19H,FA25A,QA29K,EA33R", 7, 0.747160, True],
        ["1K8R_A_B", "DA38E", 1, 3.437662, True],
        ["4G0N_A_B", "DA38A", 1, 2.864766, True],
        ["1LFD_A_B", "DA38K", 1, -0.994115, True],
        ["1GRN_A_B", "DA38E", 1, 0.389867, True],
        ["1K8R_A_B", "DA38N", 1, 1.896232, True],
        ["1LFD_A_B", "DA38A", 1, -0.442452, True],
        ["1K8R_A_B", "DA38A", 1, 1.712186, True]
    ]
    expected_data = [
        ["1K8R_A_B", "AA38E", 1, 1.725476, False],
        ["1K8R_A_B", "AA38N", 1, 0.184046, False],
        ["1K8R_A_B", "EA38A", 1, -1.725476, False],
        ["1K8R_A_B", "EA38N", 1, -1.541430, False],
        ["1K8R_A_B", "NA38A", 1, -0.184046, False],
        ["1K8R_A_B", "NA38E", 1, 1.541430, False],
        ["1LFD_A_B", "AA38K", 1, -0.551663, False],
        ["1LFD_A_B", "KA38A", 1, 0.551663, False]
    ]
    columns =  ["complex", "mutstr", "num_muts", "ddG", "original"]

    df = pd.DataFrame(data, columns=columns)
    df_exp = pd.DataFrame(expected_data, columns=columns)
    df_per = permutations(df, "complex", "mutstr", "ddG", "num_muts")

    exp_sort = df_exp.sort_values(by=["mutstr"], ignore_index=True)
    per_sort = df_per.sort_values(by=["mutstr"], ignore_index=True)

    assert per_sort.equals(exp_sort)

