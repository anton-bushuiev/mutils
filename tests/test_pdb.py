import shutil

from mutils.pdb import load_pdb, merge_chains, n_atoms, n_residues
from mutils.definitions import *


def _test_merge_chains(skempi2_df):
    tmpdir = MUTILS_SKEMPI2_DIR / 'tmp'
    tmpdir.mkdir(exist_ok=True)

    for _, row in skempi2_df.iterrows():
        # Merge
        pdbid = row['PDB Id']
        partners = [row['Partner 1'], row['Partner 2']]
        inpath = MUTILS_SKEMPI2_DIR / 'PDBs' / f'{pdbid}.pdb'
        outpath = tmpdir / f'{pdbid}.pdb'
        merge_chains(inpath, partners, outpath, verbose=False)

        # Read initial and new models from .pdb files
        model_in = load_pdb(inpath, model_id=0, verbose=False)
        model_out = load_pdb(outpath, model_id=0, verbose=False)

        # Test correct chains
        old_chains = \
            set(map(lambda chain: chain.get_id(), model_in.get_chains()))
        new_chains = \
            set(map(lambda chain: chain.get_id(), model_out.get_chains()))
        assert new_chains == \
            old_chains - set(partners[0][1:]) - set(partners[1][1:])

        # Test same number of residues and atoms in partners
        assert sum([n_residues(model_in[chain]) for chain in partners[0]]) == \
            n_residues(model_out[partners[0][0]])
        assert sum([n_atoms(model_in[chain]) for chain in partners[1]]) == \
            n_atoms(model_out[partners[1][0]])

    shutil.rmtree(tmpdir)


def test_merge_chains_random(skempi2):
    _test_merge_chains(skempi2.sample(10))


def test_merge_chains_big_multimeric(skempi2_big_multimeric_partners):
    _test_merge_chains(skempi2_big_multimeric_partners.sample(10))
