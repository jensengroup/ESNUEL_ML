### IMPORT MODULES ###
import os
import gzip
import hashlib
import argparse
import numpy as np
import lightgbm as lgb
from rdkit import Chem
from DescriptorCreator.GraphChargeShell import GraphChargeShell
### END ###

### LOAD MAA MODEL ###
lgb_model_path = 'models/elec/SMI2GCS_3_cm5_model.txt.gz'
with gzip.open(lgb_model_path, mode="rt") as file:
    lgb_model_str = file.read()
lgb_model = lgb.Booster(model_str=lgb_model_str, silent=True)
### END ###


def parse_args():
    """
    Argument parser so this can be run from the command line
    """
    parser = argparse.ArgumentParser(description='Run ESNUEL predictions from the command line')
    parser.add_argument('-s', '--smi', default='CNC(=O)Oc1cccc2c1OC(C)(C)C2',
                        help='SMILES input for ESNUEL predictions')
    return parser.parse_args()


def remove_identical_atoms(rdkit_mol, atom_list):
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(atom_list):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)
    
    atom_list = np.array(atom_list)[idx_list].tolist()
    return atom_list


def find_identical_atoms(rdkit_mol, atom_list):
    len_list = len(atom_list)
    
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(rdkit_mol.GetAtoms()):
        if atom.GetIdx() in atom_list[:len_list]:
            sym_atoms = [int(atom_idx) for atom_idx, ranking in enumerate(atom_rank) if ranking == atom_rank[idx] and atom_idx not in atom_list] 
            atom_list.extend(sym_atoms)
    return atom_list


def gen_atom_desc(smiles, atom_sites, n_shells=3):

    name = hashlib.md5(smiles.encode()).hexdigest() # SMILES MUST BE CANONICALIZED

    file_path = f'calculations/{name}/xtb.out'
    if os.path.isfile(file_path):
        rdkit_molHs = Chem.AddHs(Chem.MolFromSmiles(smiles))

        # Get CM5 charges from output and append CM5 charges to RDKit mol object 
        cm5_list = []
        natoms = int(rdkit_molHs.GetNumAtoms())
        with open(file_path, 'r') as output:
            for line_idx, line in enumerate(output):
                if 'Mulliken/CM5' in line:
                    start = line_idx + 1
                    endindex = start + natoms
                    for i in range(start, endindex):
                        line = output.readline(i)
                        cm5_atom = float(line.split()[2])
                        cm5_list.append(cm5_atom)
                    break
        
        desc_generator.rdkit_mol = rdkit_molHs
        desc_generator.cm5_list = cm5_list
    else:
        cm5_list = desc_generator.calc_CM5_charges(smiles, name=name, optimize=False, save_output=True)
    
    descriptor_vector, mapper_vector = desc_generator.create_descriptor_vector(atom_sites, n_shells=n_shells, max_neighbors=4, use_cip_sort=True)
    return descriptor_vector


def run_preds(smiles):

    smarts_of_interest = Chem.MolFromSmarts('[#6]([#6,#7])(=[OX1])[#8][#6]') # Esters, "Amides"

    rdkit_mol = Chem.MolFromSmiles(smiles)
    sites_of_interest = list(set(np.array(rdkit_mol.GetSubstructMatches(smarts_of_interest))[:,0]))
    sites_of_interest = remove_identical_atoms(Chem.AddHs(rdkit_mol), sorted(sites_of_interest))

    descriptor_vectors = gen_atom_desc(smiles, sites_of_interest, n_shells=3)
    preds = lgb_model.predict(descriptor_vectors, num_iteration=lgb_model.best_iteration)
    
    atom_sites = find_identical_atoms(rdkit_mol, [sites_of_interest[np.argmax(preds)]])
    pred_maa = float(max(preds))
    
    return atom_sites, pred_maa


# For command line use
if __name__ == "__main__":
    args = parse_args()

    desc_generator = GraphChargeShell()

    smiles = Chem.MolToSmiles(Chem.MolFromSmiles(args.smi), isomericSmiles=True) # canonicalize input smiles
    atom_sites, pred_maa = run_preds(smiles)

    print('SMILES:', smiles)
    print('Atom indices:', ','.join([str(a) for a in atom_sites]))
    print('Pred. MAA:', pred_maa)