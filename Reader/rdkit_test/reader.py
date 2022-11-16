from rdkit import Chem
from rdkit.Chem import AllChem
from models import Substance
import os

database_dir = './data/'

dir = os.fsencode(database_dir)

for file in os.listdir(dir):
    filename = os.fsdecode(file)
    if filename.endswith(".sdf"):
        substances = Chem.SDMolSupplier(database_dir + str(filename))
        for substance in substances:
            subst = Substance()
            subst.PubChem_CID = substance.GetProp('PUBCHEM_COMPOUND_CID')
            subst.PubChem_Name = substance.GetProp('PUBCHEM_IUPAC_NAME')
            subst.PubChem_Formula = substance.GetProp('PUBCHEM_MOLECULAR_FORMULA')
            subst.PubChem_Mass = substance.GetProp('PUBCHEM_MOLECULAR_WEIGHT')
            subst.save()
