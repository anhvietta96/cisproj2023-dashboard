from django.shortcuts import render
import os
from glob import glob
from pathlib import Path
from .moleculeClass import moleculeClass
from dashboard.settings import MEDIA_ROOT


def main_compound_view(request):
    properties = ['PUBCHEM_COMPOUND_CID',
                  'PUBCHEM_CACTVS_HBOND_ACCEPTOR',
                  'PUBCHEM_CACTVS_HBOND_DONOR',
                  'PUBCHEM_IUPAC_INCHI',
                  'PUBCHEM_MOLECULAR_WEIGHT',
                  'PUBCHEM_MOLECULAR_FORMULA',
                  'PUBCHEM_XLOGP3']

    mol_properties = []
    mol_names = [Path(file).stem for file in
                 glob(os.path.join(MEDIA_ROOT, '*.sdf'))]

    for mol_name in mol_names:
        mol_instance = moleculeClass("media/", mol_name)
        mol_instance.loadMolecule()
        for property_dict in mol_instance.getPropDictIterator():
            temp = [mol_name]
            for prop in properties:
                temp.append((prop, property_dict.get(prop)))
            mol_properties.append(temp)

    return render(request, 'compounds/compounds.html',
                  {'mols_with_properties': mol_properties})
