from django.shortcuts import render
import glob
from moleculeClass import moleculeClass
from pathlib import Path


def main_compound_view(request):
    properties = ['PUBCHEM_COMPOUND_CID',
                  'PUBCHEM_CACTVS_HBOND_ACCEPTOR',
                  'PUBCHEM_CACTVS_HBOND_DONOR',
                  'PUBCHEM_IUPAC_INCHI',
                  'PUBCHEM_MOLECULAR_WEIGHT',
                  'PUBCHEM_MOLECULAR_FORMULA']

    # ToDO: Use MEDIA_ROOT
    mol_names = [Path(file).stem for file in glob.glob("media/*")]
    mol_properties = []

    for mol_name in mol_names:
        mol_instance = moleculeClass("media/", mol_name)
        mol_instance.loadMolecule()
        property_dict = mol_instance.getPropDict()
        temp = [mol_name]
        for prop in properties:
            temp.append((prop, property_dict[prop]))
        mol_properties.append(temp)

    return render(request, 'compounds/compounds.html',
                  {'mols_with_properties': mol_properties})
