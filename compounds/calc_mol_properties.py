from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, inchi
from .models import Molecule
"""
In this file we iterate over a supplier and calculate all molecular properties of interest
and saves them into a dictionary
"""

class MoleculeIteratorClass:
    """
    Class to iterate over a whole supplier with molecules
    :param: filename: directory and filename of Molecule File, supplier: Supplier of Molecule File
    :return: none
    :rtype: none
    """
    def __init__(self, filename: str, supplier=Chem.SDMolSupplier):
        self.filename = filename
        self.supplier = supplier
        self.mol_list = []

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
                    continue

                mol_props = MoleculeProperties(mol)
                mol_props.save_molecule()
                self.mol_list.append(mol_props)

    def printall(self):
        """
        Dummy-Function to print out molecules
        :param: self
        :return: none
        :rtype: none
        """
        for mol in self.mol_list:
            print(mol)


class MoleculeProperties:
    """
    Class to calculate molecule properties
    :param: dir: mol: molecules of interest
    :return: none
    :rtype: none
    """
    def __init__(self, mol):
        self.mol = mol

        # primary key
        self.inchi_key = inchi.MolToInchiKey(self.mol)

        # properties for Lipinksi's rule of five
        self.molecular_weight = round(Descriptors.MolWt(self.mol), 4)
        self.log_p = round(Descriptors.MolLogP(self.mol), 6)
        self.num_h_acceptors = Lipinski.NumHAcceptors(self.mol)
        self.num_h_donors = Lipinski.NumHDonors(self.mol)

        # additional properties
        self.rotatable_bonds = Lipinski.NumRotatableBonds(self.mol)
        self.molecule_formula = rdMolDescriptors.CalcMolFormula(self.mol)

    def return_dict(self):
        """
        Returns a dictory of molecule properties
        :param: none
        :return: dictionary of molecular discriptors
        :rtype: dict
        """
        return {
            'inchi_key': self.inchi_key,
            'molecular_weight': self.molecular_weight,
            'log_p': self.log_p,
            'rotatable_bonds': self.rotatable_bonds,
            'num_h_acceptors': self.num_h_acceptors,
            'num_h_donors': self.num_h_donors,
            'molecule_formula': self.molecule_formula,
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
            num_h_acceptors=self.num_h_acceptors,
            num_h_donors=self.num_h_donors,
            log_p=self.log_p,
            molecular_mass=self.molecular_weight
        )
        m.save()

    def __str__(self):
        return str(self.return_dict())


if __name__ == '__main__':
    filename = 'Beispieldateien/Drugs/Drugs.sdf'
    x = MoleculeIteratorClass(filename)
    x.iterate_over_molecules()
    x.printall()
