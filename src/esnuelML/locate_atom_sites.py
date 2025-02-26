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

import numpy as np
from rdkit import Chem

import molecule_formats as molfmt


n_smirks_dict = {'Ether': '[OX2:1]([#6;!$(C([OX2])[#7,#8,#15,#16,F,Cl,Br,I]);!$([#6]=[#8]):2])[#6;!$(C([OX2])[#7,#8,#15,#16]);!$([#6]=[#8]):3]>>[CH3][OX3+:1]([*:2])[*:3]',
                 'Ketone': '[OX1H0:1]=[#6X3:2]([#6;!$([CX3]=[CX3;!R]):3])[#6;!$([CX3]=[CX3;!R]):4]>>[CH3][OX2H0:1][#6X3+:2]([*:3])[*:4]',
                 'Amide': '[OX1:1]=[CX3;$([CX3][#6]),$([CX3H]):2][#7X3;!R:3]>>[CH3][OX2:1][CX3:2]=[#7X3+:3]',
                 'Enolate': '[#6;$([#6]=,:[#6]-[#8-]),$([#6-]-[#6]=,:[#8]):1]~[#6:2]~[#8;$([#8-]-[#6]=,:[#6]),$([#8]=,:[#6]-[#6-]):3]>>[CH3][#6+0:1][*:2]=[#8+0:3]',
                 'Aldehyde': '[OX1:1]=[$([CX3H][#6;!$([CX3]=[CX3;!R])]),$([CX3H2]):2]>>[CH3][OX2:1][#6+:2]',
                 'Imine': '[NX2;$([N][#6]),$([NH]);!$([N][CX3]=[#7,#8,#15,#16]):1]=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6]):2]>>[CH3][NX3+:1]=[*:2]',
                 'Nitranion': '[#7X2-:1]>>[CH3][#7X3+0:1]',
                 'Carbanion': '[#6-;!$([#6X1-]#[#7,#8,#15,#16]):1]>>[CH3][#6+0:1]',
                 'Nitronate': '[#6:1]=[#7+:2](-[#8-:3])-[#8-:4]>>[CH3][#6:1][#7+:2](=[#8+0:3])-[*:4]',
                 'Ester': '[OX1:1]=[#6X3;!$([#6X3][CX3]=[CX3;!R]);$([#6X3][#6]),$([#6X3H]):2][#8X2H0:3][#6;!$(C=[O,N,S]):4]>>[CH3][OX2:1][#6X3+:2][*:3][*:4]',
                 'Carboxylic acid': '[OX1:1]=[CX3;$([R0][#6]),$([H1R0]):2][$([OX2H]),$([OX1-]):3]>>[CH3][OX2:1][CX3+:2][*:3]',
                 'Amine': '[#7+0;$([N;R;!$([#7X2]);$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]),$([NX3+0;!$([#7X3][CX3;$([CX3][#6]),$([CX3H])]=[OX1])]),$([NX4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]):1]>>[CH3][#7+:1]',
                 'Cyanoalkyl/nitrile anion': '[C:1]=[C:2]=[#7X1-:3]>>[CH3][C:1][C:2]#[#7X1+0:3]',
                 'Nitrile': '[NX1:1]#[CX2;!$(CC=C=[#7X1-]);!$(CC=C):2]>>[CH3][NX2+:1]#[*:2]',
                 'Isonitrile': '[CX1-:1]#[NX2+:2]>>[CH3][CX2+0:1]#[NX2+:2]',
                 'Phenol': '[OX2H:1][$(c(c)c),$([#6X3;R](=[#6X3;R])[#6X3;R]):2]>>[CH3][OX3+1:1][*:2]', # added due to rxn100
                 'Silyl_ether': '[#8X2H0:1][#14X4:2]([!#1:3])([!#1:4])[!#1:5]>>[CH3][#8X3H0+:1][*:2]([*:3])([*:4])[*:5]', # added due to rxn100
                 'Pyridine_like_nitrogen': '[#7X2;$([nX2](:*):*),$([#7X2;R](=[*;R])[*;R]):1]>>[CH3][#7X3+:1]', # added due to rxn100
                 'anion_with_charge_minus1': '[*-:1]>>[CH3][*+0:1]', # added to capture additional sites 
                 'double_bond': '[*;!$([!X4;!#1;!#6:1])+0:1]=[*+0:2]>>[CH3][*:1]-[*+1:2]', # added to capture additional sites
                 'double_bond_neighbouratom_with_charge_plus1': '[*;!$([!X4;!#1;!#6:1])+0:1]=[*+1:2]>>[CH3][*:1]-[*+2:2]', # added to capture additional sites
                 'triple_bond': '[*;!$([!X4;!#1;!#6:1])+0:1]#[*+0:2]>>[CH3][*:1]=[*+1:2]', # added to capture additional sites 
                 'triple_bond_neighbouratom_with_charge_plus1': '[*;!$([!X4;!#1;!#6:1])+0:1]#[*+1:2]>>[CH3][*:1]=[*+2:2]', # added to capture additional sites
                 'atom_with_lone_pair': '[!X4;!#1;!#6:1]>>[CH3][*+1:1]', # added to capture additional sites
                }

# ### BEGIN SMe (OBS! Require changes to run_rxn() in molecule_formats.py) ###
# n_smirks_dict = {'Ether': '[OX2:1]([#6;!$(C([OX2])[#7,#8,#15,#16,F,Cl,Br,I]);!$([#6]=[#8]):2])[#6;!$(C([OX2])[#7,#8,#15,#16]);!$([#6]=[#8]):3]>>[CH3][#16X2][OX3+:1]([*:2])[*:3]',
#                  'Ketone': '[OX1H0:1]=[#6X3:2]([#6;!$([CX3]=[CX3;!R]):3])[#6;!$([CX3]=[CX3;!R]):4]>>[CH3][#16X2][OX2H0:1][#6X3+:2]([*:3])[*:4]',
#                  'Amide': '[OX1:1]=[CX3;$([CX3][#6]),$([CX3H]):2][#7X3;!R:3]>>[CH3][#16X2][OX2:1][CX3:2]=[#7X3+:3]',
#                  'Enolate': '[#6;$([#6]=,:[#6]-[#8-]),$([#6-]-[#6]=,:[#8]):1]~[#6:2]~[#8;$([#8-]-[#6]=,:[#6]),$([#8]=,:[#6]-[#6-]):3]>>[CH3][#16X2][#6+0:1][*:2]=[#8+0:3]',
#                  'Aldehyde': '[OX1:1]=[$([CX3H][#6;!$([CX3]=[CX3;!R])]),$([CX3H2]):2]>>[CH3][#16X2][OX2:1][#6+:2]',
#                  'Imine': '[NX2;$([N][#6]),$([NH]);!$([N][CX3]=[#7,#8,#15,#16]):1]=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6]):2]>>[CH3][#16X2][NX3+:1]=[*:2]',
#                  'Nitranion': '[#7X2-:1]>>[CH3][#16X2][#7X3+0:1]',
#                  'Carbanion': '[#6-;!$([#6X1-]#[#7,#8,#15,#16]):1]>>[CH3][#16X2][#6+0:1]',
#                  'Nitronate': '[#6:1]=[#7+:2](-[#8-:3])-[#8-:4]>>[CH3][#16X2][#6:1][#7+:2](=[#8+0:3])-[*:4]',
#                  'Ester': '[OX1:1]=[#6X3;!$([#6X3][CX3]=[CX3;!R]);$([#6X3][#6]),$([#6X3H]):2][#8X2H0:3][#6;!$(C=[O,N,S]):4]>>[CH3][#16X2][OX2:1][#6X3+:2][*:3][*:4]',
#                  'Carboxylic acid': '[OX1:1]=[CX3;$([R0][#6]),$([H1R0]):2][$([OX2H]),$([OX1-]):3]>>[CH3][#16X2][OX2:1][CX3+:2][*:3]',
#                  'Amine': '[#7+0;$([N;R;!$([#7X2]);$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]),$([NX3+0;!$([#7X3][CX3;$([CX3][#6]),$([CX3H])]=[OX1])]),$([NX4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]):1]>>[CH3][#16X2][#7+:1]',
#                  'Cyanoalkyl/nitrile anion': '[C:1]=[C:2]=[#7X1-:3]>>[CH3][#16X2][C:1][C:2]#[#7X1+0:3]',
#                  'Nitrile': '[NX1:1]#[CX2;!$(CC=C=[#7X1-]);!$(CC=C):2]>>[CH3][#16X2][NX2+:1]#[*:2]',
#                  'Isonitrile': '[CX1-:1]#[NX2+:2]>>[CH3][#16X2][CX2+0:1]#[NX2+:2]',
#                  'Phenol': '[OX2H:1][$(c(c)c),$([#6X3;R](=[#6X3;R])[#6X3;R]):2]>>[CH3][#16X2][OX3+1:1][*:2]', # added due to rxn100
#                  'Silyl_ether': '[#8X2H0:1][#14X4:2]([!#1:3])([!#1:4])[!#1:5]>>[CH3][#16X2][#8X3H0+:1][*:2]([*:3])([*:4])[*:5]', # added due to rxn100
#                  'Pyridine_like_nitrogen': '[#7X2;$([nX2](:*):*),$([#7X2;R](=[*;R])[*;R]):1]>>[CH3][#16X2][#7X3+:1]', # added due to rxn100
#                  'anion_with_charge_minus1': '[*-:1]>>[CH3][#16X2][*+0:1]', # added to capture additional sites 
#                  'double_bond': '[*;!$([!X4;!#1;!#6:1])+0:1]=[*+0:2]>>[CH3][#16X2][*:1]-[*+1:2]', # added to capture additional sites 
#                  'double_bond_neighbouratom_with_charge_plus1': '[*;!$([!X4;!#1;!#6:1])+0:1]=[*+1:2]>>[CH3][#16X2][*:1]-[*+2:2]', # added to capture additional sites 
#                  'triple_bond': '[*;!$([!X4;!#1;!#6:1])+0:1]#[*+0:2]>>[CH3][#16X2][*:1]=[*+1:2]', # added to capture additional sites 
#                  'triple_bond_neighbouratom_with_charge_plus1': '[*;!$([!X4;!#1;!#6:1])+0:1]#[*+1:2]>>[CH3][#16X2][*:1]=[*+2:2]', # added to capture additional sites 
#                  'atom_with_lone_pair': '[!X4;!#1;!#6:1]>>[CH3][#16X2][*+1:1]', # added to capture additional sites 
#                 }
# ### END SMe ###


e_smirks_dict = {'Oxonium': '[#6:1]=[O+;!$([O]~[!#6]);!$([S]*~[#7,#8,#15,#16]):2]>>[CH3][#6:1][O+0:2]',
                 'Carbocation': '[#6+:1]>>[CH3][#6+0:1]',
                 'Ketone': '[#6X3:1](=[OX1:2])([#6;!$([CX3]=[CX3;!R]):3])[#6;!$([CX3]=[CX3;!R]):4]>>[CH3][#6X4:1](-[OX1-:2])([*:3])[*:4]',
                 'Amide': '[CX3;$([CX3][#6]),$([CX3H]):1](=[OX1:2])[#7X3;!R:3]>>[CH3][CX4:1](-[OX1-:2])[*:3]',
                 'Ester': '[#6X3;!$([#6X3][CX3]=[CX3;!R]);$([#6X3][#6]),$([#6X3H]),$([#6X3][OX2H0]):1](=[OX1:2])[#8X2H0:3][#6;!$(C=[O,N,S]):4]>>[CH3][#6X4:1](-[OX1-:2])[*:3][*:4]',
                 'Iminium': '[CX3:1]=[NX3+;!$(N([#8-])[#8-]):2]>>[CH3][CX4:1]-[NX3+0:2]',
                 'Michael acceptor': '[CX3;!R:1]=[CX3:2][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]):3]>>[CH3][CX4:1]-[CX3-:2][*:3]',
                 'Imine': '[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6]):1]=[NX2;$([N][#6]),$([NH]);!$([N][CX3]=[#7,#8,#15,#16]):2]>>[CH3][CX4:1]-[NX2-:2]',
                 'Aldehyde': '[CX3;$([CX3H][#6;!$([CX3]=[CX3;!R])]),$([CX3H2]):1]=[OX1:2]>>[CH3][CX4:1]-[OX1-:2]',
                 'Anhydride': '[CX3:1](=[OX1:2])[OX2:3][CX3:4]=[OX1:5]>>[CH3][CX4:1](-[OX1-:2])[*:3][*:4]=[*:5]', # added due to rxn100
                 'Acyl Halide': '[CX3:1](=[OX1:2])[ClX1,BrX1,IX1:3]>>[CH3][CX4:1](-[OX1-:2])[*:3]', # added due to rxn100
                 'cation_with_charge_plus1': '[*+:1]>>[CH3][*+0:1]', # added to capture additional sites 
                 'double_bond': '[*+0:1]=[*+0:2]>>[CH3][*:1]-[*-1:2]', # added to capture additional sites 
                 'double_bond_neighbouratom_with_charge_plus1': '[*+0:1]=[*+1:2]>>[CH3][*:1]-[*+0:2]', # added to capture additional sites
                 'triple_bond': '[*+0:1]#[*+0:2]>>[CH3][*:1]=[*-1:2]', # added to capture additional sites 
                 'triple_bond_neighbouratom_with_charge_plus1': '[*+0:1]#[*+1:2]>>[CH3][*:1]=[*+0:2]', # added to capture additional sites
                }

# ### BEGIN SMe (OBS! Require changes to run_rxn() in molecule_formats.py) ###
# e_smirks_dict = {'Oxonium': '[#6:1]=[O+;!$([O]~[!#6]);!$([S]*~[#7,#8,#15,#16]):2]>>[CH3][#16X2][#6:1][O+0:2]',
#                  'Carbocation': '[#6+:1]>>[CH3][#16X2][#6+0:1]',
#                  'Ketone': '[#6X3:1](=[OX1:2])([#6;!$([CX3]=[CX3;!R]):3])[#6;!$([CX3]=[CX3;!R]):4]>>[CH3][#16X2][#6X4:1](-[OX1-:2])([*:3])[*:4]',
#                  'Amide': '[CX3;$([CX3][#6]),$([CX3H]):1](=[OX1:2])[#7X3;!R:3]>>[CH3][#16X2][CX4:1](-[OX1-:2])[*:3]',
#                  'Ester': '[#6X3;!$([#6X3][CX3]=[CX3;!R]);$([#6X3][#6]),$([#6X3H]),$([#6X3][OX2H0]):1](=[OX1:2])[#8X2H0:3][#6;!$(C=[O,N,S]):4]>>[CH3][#16X2][#6X4:1](-[OX1-:2])[*:3][*:4]',
#                  'Iminium': '[CX3:1]=[NX3+;!$(N([#8-])[#8-]):2]>>[CH3][#16X2][CX4:1]-[NX3+0:2]',
#                  'Michael acceptor': '[CX3;!R:1]=[CX3:2][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]):3]>>[CH3][#16X2][CX4:1]-[CX3-:2][*:3]',
#                  'Imine': '[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6]):1]=[NX2;$([N][#6]),$([NH]);!$([N][CX3]=[#7,#8,#15,#16]):2]>>[CH3][#16X2][CX4:1]-[NX2-:2]',
#                  'Aldehyde': '[CX3;$([CX3H][#6;!$([CX3]=[CX3;!R])]),$([CX3H2]):1]=[OX1:2]>>[CH3][#16X2][CX4:1]-[OX1-:2]',
#                  'Anhydride': '[CX3:1](=[OX1:2])[OX2:3][CX3:4]=[OX1:5]>>[CH3][#16X2][CX4:1](-[OX1-:2])[*:3][*:4]=[*:5]', # added due to rxn100
#                  'Acyl Halide': '[CX3:1](=[OX1:2])[ClX1,BrX1,IX1:3]>>[CH3][#16X2][CX4:1](-[OX1-:2])[*:3]', # added due to rxn100
#                  'cation_with_charge_plus1': '[*+:1]>>[CH3][#16X2][*+0:1]', # added to capture additional sites 
#                  'double_bond': '[*+0:1]=[*+0:2]>>[CH3][#16X2][*:1]-[*-1:2]', # added to capture additional sites 
#                  'double_bond_neighbouratom_with_charge_plus1': '[*+0:1]=[*+1:2]>>[CH3][#16X2][*:1]-[*+0:2]', # added to capture additional sites 
#                  'triple_bond': '[*+0:1]#[*+0:2]>>[CH3][#16X2][*:1]=[*-1:2]', # added to capture additional sites 
#                  'triple_bond_neighbouratom_with_charge_plus1': '[*+0:1]#[*+1:2]>>[CH3][#16X2][*:1]=[*+0:2]', # added to capture additional sites 
#                 }
# ### END SMe ###


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


def remove_identical_atoms_with_associated_list(rdkit_mol, atom_list, associated_list):
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(atom_list):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)
    
    atom_list = np.array(atom_list)[idx_list].tolist()
    associated_list = np.array(associated_list)[idx_list].tolist()
    return atom_list, associated_list


def find_identical_atoms(rdkit_mol, atom_list):
    len_list = len(atom_list)
    
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(rdkit_mol.GetAtoms()):
        if atom.GetIdx() in atom_list[:len_list]:
            sym_atoms = [int(atom_idx) for atom_idx, ranking in enumerate(atom_rank) if ranking == atom_rank[idx] and atom_idx not in atom_list] 
            atom_list.extend(sym_atoms)
    return atom_list


def find_identical_atoms_with_associated_list(rdkit_mol, atom_list, associated_list):
    len_list = len(atom_list)
    
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(rdkit_mol.GetAtoms()):
        if atom.GetIdx() in atom_list[:len_list]:
            sym_atoms = [int(atom_idx) for atom_idx, ranking in enumerate(atom_rank) if ranking == atom_rank[idx] and atom_idx not in atom_list] 
            atom_list.extend(sym_atoms)
            associated_list.extend([associated_list[atom_list[:len_list].index(atom.GetIdx())]]*len(sym_atoms))
    return atom_list, associated_list


def find_electrophilic_sites(rdkit_mol):
    copy_rdkit_mol = Chem.Mol(rdkit_mol, True)
    copy_rdkit_mol = Chem.AddHs(copy_rdkit_mol)
    Chem.Kekulize(copy_rdkit_mol)

    elec_sites = []
    elec_names = []
    elec_smirks = []
    for name, smirks in e_smirks_dict.items():
        smarts = smirks.split('>>')[0]
        if copy_rdkit_mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            sites = [x[0] for x in copy_rdkit_mol.GetSubstructMatches(Chem.MolFromSmarts(smarts), uniquify=False)]
            # sites = remove_identical_atoms(copy_rdkit_mol, sites)
            for site in sites:
                if site not in elec_sites:
                    elec_sites.append(site)
                    elec_names.append(name)
                    elec_smirks.append(smirks)
    
    # Remove sites with same canonical rank
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(copy_rdkit_mol, breakTies=False))
    for idx, atom in enumerate(elec_sites):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)

    elec_sites = np.array(elec_sites)[idx_list].tolist()
    elec_names = np.array(elec_names)[idx_list].tolist()
    elec_smirks = np.array(elec_smirks)[idx_list].tolist()

    return elec_sites, elec_names, elec_smirks


def find_nucleophilic_sites(rdkit_mol):
    copy_rdkit_mol = Chem.Mol(rdkit_mol, True)
    copy_rdkit_mol = Chem.AddHs(copy_rdkit_mol)
    Chem.Kekulize(copy_rdkit_mol)

    nuc_sites = []
    nuc_names = []
    nuc_smirks = []
    for name, smirks in n_smirks_dict.items():
        smarts = smirks.split('>>')[0]
        if copy_rdkit_mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            sites = [x[0] for x in copy_rdkit_mol.GetSubstructMatches(Chem.MolFromSmarts(smarts), uniquify=False)]
            # sites = remove_identical_atoms(copy_rdkit_mol, sites)
            for site in sites:
                if site not in nuc_sites:
                    nuc_sites.append(site)
                    nuc_names.append(name)
                    nuc_smirks.append(smirks)
    
    # Remove sites with same canonical rank
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(copy_rdkit_mol, breakTies=False))
    for idx, atom in enumerate(nuc_sites):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)

    nuc_sites = np.array(nuc_sites)[idx_list].tolist()
    nuc_names = np.array(nuc_names)[idx_list].tolist()
    nuc_smirks = np.array(nuc_smirks)[idx_list].tolist()

    return nuc_sites, nuc_names, nuc_smirks


def find_nucleophilic_sites_and_generate_MCAproducts(rdkit_mol):
    nuc_prods = []
    nuc_smis = []
    nuc_sites = []
    nuc_names = []
    for name, smirks in n_smirks_dict.items():

        product_mols, product_smis, sites = molfmt.run_rxn(Chem.AddHs(rdkit_mol), smirks)
            
        for product_mol, product_smi, site in zip(product_mols, product_smis, sites):
            # if product_smi not in nuc_smis:
            if site not in nuc_sites:
                nuc_prods.append(product_mol)
                nuc_smis.append(product_smi)
                nuc_sites.append(site)
                nuc_names.append(name)
    
    # Remove sites with same canonical rank
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(nuc_sites):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)

    nuc_prods = np.array(nuc_prods)[idx_list].tolist()
    nuc_smis = np.array(nuc_smis)[idx_list].tolist()
    nuc_sites = np.array(nuc_sites)[idx_list].tolist()
    nuc_names = np.array(nuc_names)[idx_list].tolist()

    return nuc_prods, nuc_smis, nuc_sites, nuc_names


def find_electrophilic_sites_and_generate_MAAproducts(rdkit_mol):
    elec_prods = []
    elec_smis = []
    elec_sites = []
    elec_names = []
    for name, smirks in e_smirks_dict.items():
            
        product_mols, product_smis, sites = molfmt.run_rxn(Chem.AddHs(rdkit_mol), smirks)
            
        for product_mol, product_smi, site in zip(product_mols, product_smis, sites):
            # if product_smi not in elec_smis:
            if site not in elec_sites:
            # if site == 0: <-- KNUD PROJECT
                elec_prods.append(product_mol)
                elec_smis.append(product_smi)
                elec_sites.append(site)
                elec_names.append(name)
    
    # Remove sites with same canonical rank
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(elec_sites):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)

    elec_prods = np.array(elec_prods)[idx_list].tolist()
    elec_smis = np.array(elec_smis)[idx_list].tolist()
    elec_sites = np.array(elec_sites)[idx_list].tolist()
    elec_names = np.array(elec_names)[idx_list].tolist()

    return elec_prods, elec_smis, elec_sites, elec_names