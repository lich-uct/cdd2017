#!/usr/bin/env python3
### example exercise script for computational drug development
### UCT Prague, 143
import sys

from rdkit.Chem import AllChem as Chem
from rdkit import DataStructs
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

import numpy as np

### SDF reader
def sdf2mol(sdfpath):
    supplier = Chem.SDMolSupplier(sdfpath)
    for mol in supplier:
        if mol:
            yield mol

### load all the molecules in the SDF
molecules = [m for m in sdf2mol(sys.argv[1])]


### simple substructure search
query_pattern = Chem.MolFromSmarts('c1ccccc1')
filtered_molecules = [m for m in molecules if m.HasSubstructMatch(query_pattern)]


### simple similarity search
# generate molecule instance -- molecule Morgan fingerprint pairs
molecule__fp = [(m, Chem.GetMorganFingerprint(m, 2)) for m in molecules]
query_molecule = Chem.MolFromSmiles('CN(CCC1)[C@@H]1C2=CC=CN=C2') # nicotine
query_molecule_fp = Chem.GetMorganFingerprint(query_molecule, 2) # generate query fingerprint
# calculate fp similarity between the query fp and all other fps
molecule__query_similarity = [(m, DataStructs.TanimotoSimilarity(fp, query_molecule_fp))
                              for m, fp in molecule__fp]
# sort by generated similarity
sorted_m__qs = sorted(molecule__query_similarity, key=lambda x: x[1], reverse=True)
#print top ten matches
[print(Chem.MolToSmiles(m), m.GetProp('GENERIC_NAME'), sim) for m, sim in sorted_m__qs[:10]]


### simple diversity picking
# define our way of getting distances between fps under given indices
def fp_distance(i, j, fingerprints=[fp for m, fp in molecule__fp]):
    return 1 - DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])

# instantiate the RDKit diversity picker
picker = MaxMinPicker()
# run the picker with the fp_distance function, get 10 diverse pick indices
pickIndices = picker.LazyPick(fp_distance, len(molecule__fp), 10, seed=42)
# grab the Molecule instances that have those picked indices
picks = [molecule__fp[i][0] for i in pickIndices]
[print(Chem.MolToSmiles(m), m.GetProp('GENERIC_NAME')) for m in picks]

## now inverse it for similarity picking
# redefine the function
def fp_distance(i, j, fingerprints=[fp for m, fp in molecule__fp]):
    return DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
# run the picker again & print the picked molecules.
# They should be noticeably more similar than those in previous set
pickIndices = picker.LazyPick(fp_distance, len(molecule__fp), 10, seed=42)
picks = [molecule__fp[i][0] for i in pickIndices]
[print(Chem.MolToSmiles(m), m.GetProp('GENERIC_NAME')) for m in picks]
