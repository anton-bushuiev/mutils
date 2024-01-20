from pathlib import Path
import os
from collections import defaultdict
import urllib

import numpy as np
from sklearn.metrics import mean_squared_error
import Bio.PDB
from Bio import SeqIO
from Bio.PDB import PDBIO, PDBParser, parse_pdb_header
from Bio.PDB.Polypeptide import protein_letters_3to1

from mutils.misc import verbose


def path2id(path):
    return Path(path).stem.split('.', 1)[0]


@verbose
def load_pdb(path, structure_id='', model_id=None):
    if not structure_id:
        structure_id = path2id(path)
    structure = Bio.PDB.PDBParser().get_structure(structure_id, path)
    if model_id is None:
        return structure
    else:
        model = structure[model_id]
        return model


def remove_hetero_atoms(structure):
    for chain in structure:
        for res in list(chain):  # needs to be copied to detach correctly
            if res.id[0] != ' ':
                chain.detach_child(res.id)


def n_residues(entity):
    return len(list(entity.get_residues()))


def n_atoms(entity):
    return len(list(entity.get_atoms()))


# TODO: Biopython structure on input
# TODO: out to file if specified
@verbose
def pdb2fasta(path):
    fasta = ''
    with open(path, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            fasta += '>' + record.id + '\n'
            fasta += record.seq + '\n'

    return str(fasta)


@verbose
def get_sequences(pdb_path, model_id=0):
    model = load_pdb(pdb_path, model_id=model_id)
    remove_hetero_atoms(model)
    seqs = {}
    for chain in model:
        seq = []
        for residue in chain:
            seq.append(protein_letters_3to1[residue.resname])
        seqs[chain.id] = ''.join(seq)
    return seqs


@verbose
def split_pdb(path, out_path, structure_id=None):
    structure = PDBParser().get_structure(structure_id, path)
    pdb_chains = structure.get_chains()

    io = PDBIO()
    for chain in pdb_chains:
        io.set_structure(chain)
        filename = structure.get_id() + '_' + chain.get_id() + '.pdb'
        io.save(os.path.join(out_path, filename))


@verbose
def dist_interface(pdbpath, chains, dist=8, aC=True,
                   structure_id='', model_id=0):
    """
    Extract interaction interface from .pdb file.
    :param pdbpath: Path to .pdb file
    :param chains: Pair of lists of interacting chains.
        Example: (['A'], ['B', 'C']) for chain A interacting with B and C.
    :param dist: Maximum distance for residue (non-H atom) from the partnering
        chain to be considered a part of interface. [angstrom]
    :param structure_id: Structure id
    :param model_id: Model id in .pdb file
    :param aC: Whether to consider only alpha-Carbons or not
    :return: dict: chain -> list of residue ids
    """
    # Read .pdb file
    structure = load(pdbpath, structure_id)
    model = structure[model_id]

    # Define partnering chains
    chains = [[model[c] for c in chains[0]], [model[c] for c in chains[1]]]
    partners = {  # chain object -> list of partnering chain objects
        **{c: chains[1] for c in chains[0]},
        **{c: chains[0] for c in chains[1]}
    }

    # Build nearest-neighbors index for atoms
    atoms = Bio.PDB.Selection.unfold_entities(structure, 'A')
    ns = Bio.PDB.NeighborSearch(atoms)

    # Check if an atom is close to some partnering chain
    interface = defaultdict(list)
    for atom in atoms:
        if aC and atom.get_name() != 'CA':
            continue
        if atom.get_name() == 'H':
            continue
        neigh_chains = ns.search(atom.get_coord(), dist, level='C')
        residue = atom.get_parent()
        chain = residue.get_parent()
        for neigh_chain in neigh_chains:
            if neigh_chain in partners[chain]:
                interface[chain.get_id()].append(residue.get_id()[1])

    return dict(interface)


@verbose
def merge_chains(pdbpath, chain_groups, outfile, structure_id='', model_id=0):
    """
    Merge groups of chains into single ones. The id of a merged chain
    corresponds to the first one in the list.
    :param chain_groups: Example: [[A, B], [C, D]]
    :param outfile: Path to output .pdb file with merged chains.
        Example: outfile will contain the structure from `pdbpath` but with
        chain A comprising former A and B and C comprising C and D.
    """
    # Read .pdb file
    structure = load_pdb(pdbpath, structure_id)
    model = structure[model_id]

    for chains in chain_groups:
        # Insert residues from other chains
        last_res_id = int(model[chains[0]].get_unpacked_list()[-1].id[1])
        for chain in model:
            if chain.get_id() in chains[1:]:
                for res in chain:
                    # Modify residue id to match number of destination chain
                    res.detach_parent()
                    res_id = list(res.get_id())
                    res_id[1] = last_res_id + 1
                    last_res_id += 1
                    res.id = tuple(res_id)
                    # Add residue
                    model[chains[0]].add(res)

        # Delete other chains
        for chain in chains[1:]:
            model.detach_child(chain)

    # Save merged
    pdb_io = Bio.PDB.PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(str(outfile))


def coord_diff(s1, s2, substructure=None):
    """
    Compare coordinates of two structures of the same protein (e.g, after
    force field optimization).
    :param s1: Bio.PDB.Model
    :param s2: Bio.PDB.Model
    :param substructure: Dict (chain id -> list of residue ids) to specify
        residues to consider. If not specified, all residues are compared
    :return: Number of moved residues, number of moved atoms,
        mean deviation of moved atoms, RMSD of moved atoms
    """
    n_residues = 0
    n_residues_diff = 0
    n_atoms = 0
    n_atoms_diff = 0
    res_diff = False
    changed_coords_1, changed_coords_2 = [], []
    for chain1, chain2 in zip(s1, s2):
        for res1, res2 in zip(chain1, chain2):
            # Skip residue if it is not in specified substructure
            if substructure is not None and \
                    (chain1.id not in substructure
                     or res1.id[1] not in substructure[chain1.id]):
                continue
            n_residues += 1

            # Compare atom coordinates
            for atom1, atom2 in zip(res1, res2):  # Does not consider new atoms
                n_atoms += 1
                coords1 = list(atom1.get_vector())
                coords2 = list(atom2.get_vector())
                if coords1 != coords2:
                    n_atoms_diff += 1
                    if not res_diff:
                        res_diff = True
                        n_residues_diff += 1
                    changed_coords_1.append(coords1)
                    changed_coords_2.append(coords2)
            res_diff = False

    # Calculate ratios of moved atoms
    n_residues_diff_ratio = n_residues_diff / n_residues
    n_atoms_diff_ratio = n_atoms_diff / n_atoms

    if len(changed_coords_1) != 0:
        changed_coords_1 = np.array(changed_coords_1)
        changed_coords_2 = np.array(changed_coords_2)

        # Calculate mean deviation of moved atoms
        md = np.linalg.norm(changed_coords_1 - changed_coords_2, axis=1).mean()

        # Calculate RMSD of moved atoms
        rmsd = mean_squared_error(
            changed_coords_1.flatten(),
            changed_coords_2.flatten(),
            squared=False
        )
    else:
        md, rmsd = 0, 0

    return n_residues_diff, n_residues_diff_ratio, \
           n_atoms_diff, n_atoms_diff_ratio, \
           md, rmsd


def download_pdb(pdb_id, dir='.', path=None):
    if path is None:
        path = Path(dir) / f'{pdb_id}.pdb'
    path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(
        f'http://files.rcsb.org/download/{pdb_id}.pdb',
        str(path)
    )


def get_resolution(pdb_id):
    path = f'tmp-{pdb_id}.pdb'
    download_pdb(pdb_id, path=path)
    resolution = parse_pdb_header(path)['resolution']
    os.remove(path)
    return resolution
