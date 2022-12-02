import os.path
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, inchi, Crippen
from .models import Molecule

"""
In this file we iterate over a supplier and calculate all molecular properties 
of interest and save them into a dictionary
"""


class MoleculeIterator:
    """
    Class to iterate over all molecules in an .sdf or .smi file
    :param: filename: directory and filename of Molecule File
    :return: none
    :rtype: none
    """

    def __init__(self, filename: str, save_to_list: bool = True,
                 multithreaded: bool = False):

        self.filename = filename
        self.save_mol_list = save_to_list
        self.multithreaded = multithreaded
        self.supplier = None
        self.mol_list = []

        if not os.path.isfile(filename):
            raise ValueError

        if filename.endswith('.sdf'):
            self.supplier = Chem.SDMolSupplier if not multithreaded \
                else Chem.MultithreadedSDMolSupplier
        elif filename.endswith('.smi'):
            self.supplier = Chem.SmilesMolSupplier if not multithreaded \
                else Chem.MultithreadedSmilesMolSupplier
        else:
            raise ValueError

    def iterate_over_molecules(self):
        """
        This function iterates over all molecules of a supplier
        :param: self
        :return: none
        :rtype: none
        """
        with self.supplier(self.filename) as suppl:
            for mol in suppl:
                if not mol:
                    print("Cannot parse molecule", file=sys.stderr)
                    continue

                mol_props = MoleculeProperties(mol)
                mol_props.save_molecule()

                if self.save_mol_list:
                    self.mol_list.append(mol_props)

    def printall(self):
        """
        Dummy-Function to print out molecules
        :param: self
        :return: none
        :rtype: none
        """
        if not self.save_mol_list:
            print("Saving molecules is disabled", file=sys.stderr)

        for mol in self.mol_list:
            print(mol)


class MoleculeProperties:
    """
    Class to calculate molecule properties
    :param: mol: molecule of interest
    :return: none
    :rtype: none
    """

    def __init__(self, mol):
        self.mol = mol

        # primary key
        self.inchi_key = inchi.MolToInchiKey(self.mol)

        # properties for Lipinksi's rule of five
        self.log_p = round(Crippen.MolLogP(self.mol), 6)
        self.molecular_weight = round(Descriptors.MolWt(self.mol), 4)
        self.num_h_acceptors = Lipinski.NumHAcceptors(self.mol)
        self.num_h_donors = Lipinski.NumHDonors(self.mol)

        # additional properties
        self.molecular_formula = rdMolDescriptors.CalcMolFormula(self.mol)
        self.num_rings = rdMolDescriptors.CalcNumRings(self.mol)
        self.rotatable_bonds = Lipinski.NumRotatableBonds(self.mol)

    def return_dict(self):
        """
        Returns a dictionary of molecule properties
        :param: none
        :return: dictionary of molecular descriptors
        :rtype: dict
        """
        return {
            'inchi_key': self.inchi_key,
            'log_p': self.log_p,
            'molecular_formula': self.molecular_formula,
            'molecular_weight': self.molecular_weight,
            'num_h_acceptors': self.num_h_acceptors,
            'num_h_donors': self.num_h_donors,
            'num_rings': self.num_rings,
            'num_rotatable_bonds': self.rotatable_bonds
        }

    def save_molecule(self):
        """
        Saves molecular properties
        :param: none
        :return: none
        :rtype: none
        """
        m = Molecule(
            inchi_key=self.inchi_key,
            log_p=self.log_p,
            molecular_formula=self.molecular_formula,
            molecular_weight=self.molecular_weight,
            num_h_acceptors=self.num_h_acceptors,
            num_h_donors=self.num_h_donors,
            num_rings=self.num_rings,
            num_rotatable_bonds=self.rotatable_bonds
        )
        m.save()

    def __str__(self):
        return str(self.return_dict())


if __name__ == '__main__':
    filename = 'Beispieldateien/Drugs/Drugs.sdf'
    x = MoleculeIterator(filename)
    x.iterate_over_molecules()
    x.printall()
