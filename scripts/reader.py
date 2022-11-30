from rdkit import Chem
from compounds.models import Molecule
import os
import re
from pathlib import Path

def run():
    database_dir = './data/'
    dir = os.fsencode(database_dir)
    file_list = [f for f in Path(database_dir).glob('**/*') if f.is_file()]
    for file in file_list:
        filename = os.fsdecode(file)
        if filename.endswith(".sdf"):
            with Chem.SDMolSupplier(filename) as suppl:
                for substance in suppl:
                    if substance is None:
                        continue
                    subst = Molecule()
                    subst.inchi_key = substance.GetProp('PUBCHEM_IUPAC_INCHIKEY')
                    subst.log_p = substance.GetProp('PUBCHEM_XLOGP3')
                    subst.num_h_acceptors = substance.GetProp('PUBCHEM_CACTVS_HBOND_ACCEPTOR')
                    subst.num_h_donors = substance.GetProp('PUBCHEM_CACTVS_HBOND_DONOR')
                    subst.molecular_mass = substance.GetProp('PUBCHEM_MOLECULAR_WEIGHT')
                    subst.save()
