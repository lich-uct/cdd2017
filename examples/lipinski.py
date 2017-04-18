#!/usr/bin/env python3
### example exercise script for computational drug development
### UCT Prague, 143
import sys

import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors as rdescriptors
import numpy as np

### just wrapping the basic RDKit functionality, in case it needs to be changed later
def num_hydrogen_bond_acceptors(mol):
    return rdescriptors.CalcNumLipinskiHBA(mol)

def num_hydrogen_bond_donors(mol):
    return rdescriptors.CalcNumLipinskiHBD(mol)

def MW(mol):
    return Descriptors.MolWt(mol)

def logP(mol):
    return Descriptors.MolLogP(mol)

def TPSA(mol):
    return Descriptors.TPSA(mol)

def num_rotatable_bonds(mol):
    return Descriptors.NumRotatableBonds(mol)

def num_heavy_atoms(mol):
    return mol.GetNumHeavyAtoms()

### Lipinski checks, return True/False
def lipinski_HBA(mol):
    return num_hydrogen_bond_acceptors(mol) <= 10

def lipinski_HBD(mol):
    return num_hydrogen_bond_acceptors(mol) <= 5

def lipinski_MW(mol):
    return MW(mol) < 500

def lipinski_logP(mol):
    return logP(mol) <= 5

def lipinski_TPSA(mol):
    return TPSA(mol) < 140

def lipinski_rotatable_bonds(mol):
    return num_rotatable_bonds(mol) < 10

def lipinski_heavy_atoms(mol):
    return 20 < num_heavy_atoms(mol) < 70

### Lipinski rule violation counter
def num_lipinski_violations(mol, *args, ruleset='basic'):
    rulesets = {'basic': (lipinski_HBA, lipinski_HBD, lipinski_logP, lipinski_MW),
                'extended': (lipinski_HBA, lipinski_HBD, lipinski_logP, lipinski_MW, lipinski_TPSA, lipinski_heavy_atoms, lipinski_rotatable_bonds)}
    rulescount = len(rulesets[ruleset])
    # count of all rules - count of passed rules = count of failed rules
    return rulescount - sum([f(mol) for f in rulesets[ruleset]])

### SDF reader
def sdf2mol(sdfpath):
    supplier = Chem.SDMolSupplier(sdfpath)
    for mol in supplier:
        if mol:
            yield mol

if __name__ == "__main__":
    violations = [num_lipinski_violations(m, ruleset='extended') for m in sdf2mol(sys.argv[1])]
    print(float(sum(violations)) / len(violations))

    n, bins, patches = plt.hist(violations, bins=8, range=(0,8))
    plt.show()

