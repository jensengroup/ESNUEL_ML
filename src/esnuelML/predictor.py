# MIT License
#
# Copyright (c) 2024 Nicolai Ree
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

### IMPORT MODULES ###
import os
import sys
import gzip
import hashlib
import argparse
import numpy as np
import pandas as pd
from joblib import load

from rdkit import Chem
import lightgbm as lgb

from DescriptorCreator.GraphChargeShell import GraphChargeShell
desc_generator = GraphChargeShell()

from locate_atom_sites import find_nucleophilic_sites, find_electrophilic_sites
from molecule_drawer import generate_structure, generate_output_tables, html_output

base_dir = os.path.dirname(os.path.realpath(__file__)).replace('/src/esnuelML', '')
### END ###


### LOAD MCA and MAA MODELS ###
# Property prediction models:
nuc_model_path = os.path.join(base_dir, 'src/esnuelML/models/nuc/SMI2GCS_3_cm5_model.txt.gz')
with gzip.open(nuc_model_path, mode="rt") as file:
    nuc_model_str = file.read()
nuc_model = lgb.Booster(model_str=nuc_model_str, silent=True)

elec_model_path = os.path.join(base_dir, 'src/esnuelML/models/elec/SMI2GCS_3_cm5_model.txt.gz')
with gzip.open(elec_model_path, mode="rt") as file:
    elec_model_str = file.read()
elec_model = lgb.Booster(model_str=elec_model_str, silent=True)

# Error estimation models:
nuc_RFmodel = load(os.path.join(base_dir, 'src/esnuelML/models/nuc/RF_SMI2GCS_3_cm5.joblib'))
elec_RFmodel = load(os.path.join(base_dir, 'src/esnuelML/models/elec/RF_SMI2GCS_3_cm5.joblib'))
### END ###


def parse_args():
    """
    Argument parser so this can be run from the command line
    """
    parser = argparse.ArgumentParser(description='Run esnuelML predictions from the command line')
    parser.add_argument('-s', '--smi', default='C[C+:20](C)CC(C)(C)C1=C(C=CC(=C1)Br)[OH:10]',
                        help='SMILES input for esnuelML predictions')
    return parser.parse_args()


def estimate_uncertainty(RFmodel, descs):
    tree_predictions = []
    for tree in RFmodel.estimators_:
        tree_predictions.append(tree.predict(descs))
    return np.std(tree_predictions, axis=0)


def run_MAA_and_MCA_predictions(name, smiles):

    #smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) # SMILES should already be canonicalized
    
    # Calculate CM5 atomic charges
    cm5_list = desc_generator.calc_CM5_charges(smiles, name=name, optimize=False, save_output=True)

    # Save structures in SDF format
    writer = Chem.rdmolfiles.SDWriter(desc_generator.xyz_file_path.replace('.xyz', '.sdf'))
    writer.write(desc_generator.rdkit_mol)
    writer.close()

    # Locate MAA and MCA sites
    elec_sites, elec_names, elec_smirks = find_electrophilic_sites(desc_generator.rdkit_mol) # elec_sites: site in reactant, elec_names: name of functional group
    nuc_sites, nuc_names, nuc_smirks = find_nucleophilic_sites(desc_generator.rdkit_mol) # nuc_sites: site in reactant, nuc_names: name of functional group

    # Generate atomic descriptors
    common_sites = list(set(elec_sites) & set(nuc_sites))
    common_descriptor_vectors, _ = desc_generator.create_descriptor_vector(common_sites, n_shells=3, max_neighbors=4, use_cip_sort=True)

    elec_sites_without_common = list(set(elec_sites) - set(common_sites))
    elec_without_common_descriptor_vectors, _ = desc_generator.create_descriptor_vector(elec_sites_without_common, n_shells=3, max_neighbors=4, use_cip_sort=True)
    
    nuc_sites_without_common = list(set(nuc_sites) - set(common_sites))
    nuc_without_common_descriptor_vectors, _ = desc_generator.create_descriptor_vector(nuc_sites_without_common, n_shells=3, max_neighbors=4, use_cip_sort=True)

    # Predict MAA values
    elec_descriptor_vectors = []
    for site in elec_sites:
        if site in common_sites:
            idx = common_sites.index(site)
            elec_descriptor_vectors.append(common_descriptor_vectors[idx])
        else:
            idx = elec_sites_without_common.index(site)
            elec_descriptor_vectors.append(elec_without_common_descriptor_vectors[idx])

    if len(elec_descriptor_vectors):
        elec_preds = elec_model.predict(elec_descriptor_vectors, num_iteration=elec_model.best_iteration)
        elec_errors = 0.6028123658016926 * estimate_uncertainty(elec_RFmodel, elec_descriptor_vectors) + 4.044501905413345 # If errors are above cutoff of 30.0 kJ/mol, we recommend running QM calculations
    else:
        elec_preds = []
        elec_errors = []

    # Predict MCA values
    nuc_descriptor_vectors = []
    for site in nuc_sites:
        if site in common_sites:
            idx = common_sites.index(site)
            nuc_descriptor_vectors.append(common_descriptor_vectors[idx])
        else:
            idx = nuc_sites_without_common.index(site)
            nuc_descriptor_vectors.append(nuc_without_common_descriptor_vectors[idx])

    if len(nuc_descriptor_vectors):
        nuc_preds = nuc_model.predict(nuc_descriptor_vectors, num_iteration=nuc_model.best_iteration)
        nuc_errors = 0.6140308723543682 * estimate_uncertainty(nuc_RFmodel, nuc_descriptor_vectors) + 4.3372205727670305 # If errors are above cutoff of 25.0 kJ/mol, we recommend running QM calculations
    else:
        nuc_preds = []
        nuc_errors = []
    
    return elec_sites, elec_names, elec_preds, elec_errors, nuc_sites, nuc_names, nuc_preds, nuc_errors


def pred_MAA_and_MCA(reac_smis: str, name: str):

    reac_smis = reac_smis.split('.')
    reac_mols = [Chem.MolFromSmiles(smi) for smi in reac_smis]
    reac_smis = [Chem.MolToSmiles(reac_mol) for reac_mol in reac_mols]

    # Locate MAA and MCA sites and predict MAA and MCA values in kJ/mol
    elec_sites_list = []
    elec_names_list = []
    MAA_values = []
    MAA_estimated_error = []
    elec_sdfpath_structures = []

    nuc_sites_list = []
    nuc_names_list = []
    MCA_values = []
    MCA_estimated_error = []
    nuc_sdfpath_structures = []

    for i, smiles in enumerate(reac_smis):

        elec_sites, elec_names, elec_preds, elec_errors, nuc_sites, nuc_names, nuc_preds, nuc_errors = run_MAA_and_MCA_predictions(f'{name}_{i}', smiles)

        elec_sites_list.append(elec_sites)
        elec_names_list.append(elec_names)
        MAA_values.append(elec_preds)
        MAA_estimated_error.append(elec_errors)
        elec_sdfpath_structures.append([f'{name}_{i}.sdf']*len(elec_sites))
        
        nuc_sites_list.append(nuc_sites)
        nuc_names_list.append(nuc_names)
        MCA_values.append(nuc_preds)
        MCA_estimated_error.append(nuc_errors)
        nuc_sdfpath_structures.append([f'{name}_{i}.sdf']*len(nuc_sites))
    
    ### Draw the output ###
    result_svg = generate_structure(reac_mols, elec_sites_list, MAA_values, nuc_sites_list, MCA_values, molsPerRow=2)

    df_elec = generate_output_tables(reac_mols, elec_names_list, MAA_values, MAA_estimated_error, elec_sites_list, elec_sdfpath_structures, MAA_or_MCA='MAA', QM_or_ML='ML')
    df_elec = df_elec.rename(columns={'Error Log (Reactant, Product)': 'Reactant'})
    df_nuc = generate_output_tables(reac_mols, nuc_names_list, MCA_values, MCA_estimated_error, nuc_sites_list, nuc_sdfpath_structures, MAA_or_MCA='MCA', QM_or_ML='ML')
    df_nuc = df_nuc.rename(columns={'Error Log (Reactant, Product)': 'Reactant'})

    result_output = html_output('.'.join(reac_smis), result_svg, df_elec, df_nuc)
    
    fd = open(f'{base_dir}/desc_calcs/{name}/{name}.html','w')
    fd.write(result_output)
    fd.close()
    ### END ###

    return df_elec, df_nuc, elec_names_list, MAA_values, elec_sites_list, nuc_names_list, MCA_values, nuc_sites_list


if __name__ == "__main__":

    args = parse_args() # Obtain CLI

    smiles = Chem.MolToSmiles(Chem.MolFromSmiles(args.smi), isomericSmiles=True) # canonicalize input smiles
    name = hashlib.md5(smiles.encode()).hexdigest() # SMILES MUST BE CANONICALIZED

    df_elec, df_nuc, elec_names_list, MAA_values, elec_sites_list, nuc_names_list, MCA_values, nuc_sites_list = pred_MAA_and_MCA(smiles, name)

    print('Name:', name)
    print('SMILES:', smiles)
    print(df_elec.drop(columns=['Reactant']))
    print(df_nuc.drop(columns=['Reactant']))
