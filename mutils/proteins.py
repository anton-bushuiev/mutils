from io import StringIO

import numpy as np
import pandas as pd
from Bio import SeqIO


# Based on https://en.wikipedia.org/wiki/Amino_acid
AMINO_ACID_NAMES = [
    'Alanine', 'Arginine', 'Asparagine', 'Aspartate', 'Cysteine',
    'Glutamine', 'Glutamate', 'Glycine', 'Histidine', 'Neutral',
    'Isoleucine', 'Leucine', 'Lysine', 'Methionine', 'Phenylalanine',
    'Proline', 'Serine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valin'
]
AMINO_ACID_CODES_3 = [
    'Ala', 'Arg', 'Asn', 'Asp', 'Cys',
    'Gln', 'Glu', 'Gly', 'His', 'Ile',
    'Leu', 'Lys', 'Met', 'Phe', 'Pro',
    'Ser', 'Thr', 'Trp', 'Tyr', 'Val'
]
AMINO_ACID_CODES_3_UPPER = list(map(
    lambda code: code.upper(), AMINO_ACID_CODES_3
))
AMINO_ACID_CODES_1 = [
    'A', 'R', 'N', 'D', 'C',
    'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P',
    'S', 'T', 'W', 'Y', 'V'
]
MISSING_AMINO_ACID_CODE = '-'


def parse_fasta(fasta: str) -> dict:
    """
    :param fasta: string in FASTA format
    :return: mapping of chain ids to amino acid sequences
    """
    fasta_sequences = SeqIO.parse(StringIO(fasta), 'fasta')
    sequences = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        name = name.split(':')[-1]
        sequences[name] = sequence

    return sequences


def bio_pssm_to_df(pssm, alphabet=AMINO_ACID_CODES_1):
    df = []
    for row in pssm.pssm:
        value_dict = row[1]
        df_row = {}
        for letter in alphabet:
            if letter in value_dict:
                df_row[letter] = value_dict[letter]
            else:
                df_row[letter] = np.nan
        df.append(df_row)
    return pd.DataFrame(df)


def calc_dG(affinity, temperature=None):
    """
    From https://life.bsc.es/pid/skempi2/info/faq_and_help#5
    """
    if temperature is None:
        temperature = 273.15 + 25.0
    return (8.314/4184) * temperature * np.log(affinity)
