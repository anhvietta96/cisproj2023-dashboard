import os.path
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, inchi, Crippen
from .models import Molecule, MoleculeSet

"""
In this file we iterate over a supplier and calculate all molecular properties 
of interest
"""


def insert_into_new_set(dirname, set_name='') -> None:
    """
    Adds the compounds from MEDIA_ROOT into the DB and links them to a new
    MoleculeSet with name set_name
    :param dirname: name of the directory which includes SDF or SMI files
    :param set_name: name of the new MoleculeSet
    :return:
    """
    file_iterator = FileIterator(dirname)
    file_iterator.iterate_over_files()

    mol_set = MoleculeSet(set_name=set_name)
    mol_set.save()

    mol_set.molecules.add(*file_iterator.get_mol_list())
    mol_set.save()


class FileIterator:
    """
    Class for the iteration over all files in a directory,
    adding the molecules into the DB.
    """

    def __init__(self, dirname: str):
        self.dirname = dirname
        self.mol_list = []

        if not os.path.isdir(dirname):
            print(f"{dirname} is not a valid directory", file=sys.stderr)
            raise ValueError

    def iterate_over_files(self):
        for root, dirs, files in os.walk(self.dirname):
            for filename in files:
                path_to_file = os.path.join(root, filename)

                try:
                    mol_iterator = MoleculeIterator(path_to_file)
                except ValueError:
                    print(f"Skipped {path_to_file}", file=sys.stderr)
                    continue

                mol_iterator.iterate_over_molecules()
                self.mol_list.extend(mol_iterator.get_mol_list())

    def get_mol_list(self) -> list[Molecule]:
        """
        Returns a list of all Molecule instances from self.dirname
        :return: list
        """
        return self.mol_list

    def add_mol_list_to_set(self, molecule_set: MoleculeSet) -> None:
        """
        Insert all molecules into the given MoleculeSet
        :return: None
        """
        molecule_set.molecules.add(*self.mol_list)
        molecule_set.save()


class MoleculeIterator:
    """
    Class to iterate over all molecules in an .sdf or .smi file
    :param: filename: path to file (incl. filename) of Molecule File
    :return: none
    :rtype: none
    """

    def __init__(self, filename: str,
                 save_model_instances: bool = True,
                 print_mol_dicts: bool = False,
                 multithreaded: bool = False):

        self.filename = filename
        self.supplier = None
        self.multithreaded = multithreaded
        self.print_mol_dicts = print_mol_dicts
        self.save_model_list = save_model_instances
        self.mol_list = []

        if not os.path.isfile(filename):
            print(f"{filename} is not a valid file", file=sys.stderr)
            raise ValueError

        if filename.endswith('.sdf'):
            self.supplier = Chem.SDMolSupplier if not multithreaded \
                else Chem.MultithreadedSDMolSupplier
        elif filename.endswith('.smi'):
            self.supplier = Chem.SmilesMolSupplier if not multithreaded \
                else Chem.MultithreadedSmilesMolSupplier
        else:
            print(f"{filename} is not an SDF or SMI file! ", file=sys.stderr)
            raise ValueError

    def iterate_over_molecules(self) -> None:
        """
        This function iterates over all molecules in the provided filename
        and saves the molecules into the DB
        :param: self
        :return: none
        :rtype: none
        """
        with self.supplier(self.filename) as suppl:
            for mol in suppl:
                if not mol:
                    print("Cannot parse molecule", file=sys.stderr)
                    continue

                mol_props_instance = MoleculeProperties(mol)
                mol_instance = mol_props_instance.save_molecule()

                if self.save_model_list:
                    self.mol_list.append(mol_instance)
                if self.print_mol_dicts:
                    print(mol_instance)

    def get_mol_list(self) -> list[Molecule]:
        """
        Returns a list of all Molecule instances from self.filename
        :return: list
        """
        return self.mol_list


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

    def save_molecule(self) -> Molecule:
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
        return m

    def __str__(self):
        return str(self.return_dict())


if __name__ == '__main__':
    drugs_filepath = 'Beispieldateien/Drugs/Drugs.sdf'
    x = MoleculeIterator(drugs_filepath)
    x.iterate_over_molecules()
