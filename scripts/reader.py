from rdkit import Chem
from rdkit.Chem import AllChem
from Reader.models import Substance
import os

'''os.environ.setdefault('DJANGO_SETTINGS_MODULE','testsite.settings')
django.setup()'''

def run():
    database_dir = './scripts/data/'
    dir = os.fsencode(database_dir)
    for file in os.listdir(dir):
        filename = os.fsdecode(file)
        if filename.endswith(".sdf"):
            with Chem.SDMolSupplier(database_dir + filename) as suppl:
                for substance in suppl:
                    if substance is None: continue
                    subst = Substance()
                    subst.PubChem_CID = substance.GetProp('PUBCHEM_COMPOUND_CID')
                    subst.PubChem_Name = substance.GetProp('PUBCHEM_IUPAC_NAME')
                    subst.PubChem_Alias = substance.GetProp('PUBCHEM_IUPAC_NAME')
                    subst.PubChem_Formula = substance.GetProp('PUBCHEM_MOLECULAR_FORMULA')
                    subst.PubChem_Mass = substance.GetProp('PUBCHEM_MOLECULAR_WEIGHT')
                    subst.save()
