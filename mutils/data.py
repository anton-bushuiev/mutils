import os
from collections import defaultdict

import numpy as np
import pandas as pd
import requests

from mutils.mutations import Mutation
from mutils.proteins import calc_dG
from mutils.definitions import MUTILS_SKEMPI2_DIR


def load_SKEMPI2(
    table_path=MUTILS_SKEMPI2_DIR / 'skempi_v2.csv',
    pdbs_path=MUTILS_SKEMPI2_DIR / 'PDBs',
    preprocess=True,
    dG=True,
    ddG=True,
    use_temperature_col=False,
    pdb_id=True,
    partners=True,
    drop_no_ddG=True,
    drop_nmr=True,
    single_point_only=False,
    mutations_as_list=False,
    merged_chains=False,
    raw=False
):
    # TODO Download SKEMPI2 if files not found

    # Read table
    df = pd.read_csv(table_path, sep=';', low_memory=False)

    if raw:
        return df

    # Preprocess
    if preprocess:
        # Convert types
        df['#Pdb'] = df['#Pdb'].astype(str)

        # Fill missing values
        df['Hold_out_type'] = df['Hold_out_type'].fillna('Other')

    # Calculate binding energies and binding energy changes
    if ddG or dG:
        df['Affinity_mut (M)'] = pd.to_numeric(df['Affinity_mut_parsed'],
                                               errors='coerce')
        df['Affinity_wt (M)'] = pd.to_numeric(df['Affinity_wt_parsed'],
                                              errors='coerce')
        temperature = df['Temperature'].str[:3].astype(float) if use_temperature_col else None
        dg_mut = calc_dG(df['Affinity_mut (M)'], temperature=temperature)
        dg_wt = calc_dG(df['Affinity_wt (M)'], temperature=temperature)
        if dG:
            df['dG_mut'] = dg_mut
            df['dG_wt'] = dg_wt
        if ddG:
            df['ddG'] =  dg_mut - dg_wt

    # Drop NMR
    if drop_nmr:
        df = df[df['#Pdb'] != '1KBH_A_B']

    # Drop rows containing no ddG
    if drop_no_ddG:
        df = df[df['ddG'].notna()]

    # Drop multi-point
    if single_point_only:
        df = df[df['Mutation(s)_cleaned'].str.count(',') == 0]

    # Parse `#Pdb` column
    if pdb_id or partners or merged_chains:
        parsed_pdb = df['#Pdb'].str.split('_').tolist()
        parsed_pdb = np.array(parsed_pdb)

        # Extract PDB id
        if pdb_id:
            df['PDB Id'] = parsed_pdb[:, 0]

        # Extract partners
        if partners or merged_chains:
            df['Partner 1'] = parsed_pdb[:, 1]
            df['Partner 2'] = parsed_pdb[:, 2]

    # Rename merged chains and add columns with merged partners' ids
    if merged_chains:
        def rename_mut_chains_after_merge(row):
            mut = row['Mutation(s)_cleaned']
            partner1 = row['Partner 1']
            partner2 = row['Partner 2']
            rename_dict = {
                **{p: partner1[0] for p in partner1},
                **{p: partner2[0] for p in partner2}
            }
            return str(Mutation(mut).rename_chains(rename_dict))

        df['Mutation(s)_cleaned_merged'] = \
            df.apply(rename_mut_chains_after_merge, axis=1)
        if partners:
            df['Partner 1 merged'] = df['Partner 1'].str[0]
            df['Partner 2 merged'] = df['Partner 2'].str[0]
        else:
            df = df.drop(columns=['Partner 1', 'Partner 2'])

    # Split mutations to lists
    if mutations_as_list:
        mutlists = df['Mutation(s)_cleaned'].str.replace(' ', '').str.split(',')
        df['Mutation(s)_cleaned'] = mutlists

    # TODO Read structures
    pdbs = None

    return df, pdbs


class SCOP:
    def __init__(self):
        self.df = pd.read_csv(
            'http://scop.mrc-lmb.cam.ac.uk/files/scop-cla-latest.txt',
            sep=' ',
            skiprows=5
        )

    def ancestry(self, pdbid=None, scopid=None):
        """
        :return: dict: SCOP node related to specified id -> ancestry
        """
        df = self.df

        # Get SCOP ids
        if pdbid is not None:
            scopids = df[df['FA-DOMID'] == pdbid]['#'].tolist()
        else:
            scopids = [scopid]

        # Query
        ancestry = defaultdict(list)
        for scopid in scopids:
            url = f'https://scop.mrc-lmb.cam.ac.uk/api/ancestry/{scopid}'
            response = requests.get(url)
            if response.ok:
                for node, val in response.json()['lineage']['nodes'].items():
                    ancestry[scopid].append(val)

        return dict(ancestry)

    @staticmethod
    def match(ancestry_a, ancestry_b):
        """
        :param ancestry_a:
        :param ancestry_b:
        :return: set of match types (e.g. match in superfamily, fold, ...)
        """
        classes_a = defaultdict(set)
        for node, ancestry in ancestry_a.items():
            for entry in ancestry:
                classes_a[entry['name']].add(entry['id'])

        matches = set()
        for node, ancestry in ancestry_b.items():
            for entry in ancestry:
                if entry['id'] in classes_a[entry['name']]:
                    matches.add(entry['type'])

        return matches
