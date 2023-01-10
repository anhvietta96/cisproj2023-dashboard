"""
File handling classes and functions for SDF files
"""
import os.path
import sys
from rdkit import Chem
from .models import Molecule, MoleculeSet
from .calc_mol_properties import MoleculeProperties
from django.conf import settings
from typing import Optional


def insert_into_new_set(dirname: str, set_name: str) -> None:
    """
    Adds the compounds from dirname into the DB and links them to a new
    MoleculeSet with name set_name
    :param dirname: name of the directory which includes SDF files
    :param set_name: name of the new MoleculeSet
    """
    file_iterator = FileIterator(dirname)
    file_iterator.iterate_over_files()

    mol_set = MoleculeSet(set_name=set_name)
    mol_set.save()

    mol_set.molecules.add(*file_iterator.get_mol_list())
    mol_set.save()


def insert_into_existing_set(dirname: str, set_name: str) -> None:
    """
    Adds the compounds from dirname into the DB and links them to an already
    existing MoleculeSet with id set_id
    :param dirname: name of the directory which includes SDF files
    :param set_name: name of the new MoleculeSet
    """
    file_iterator = FileIterator(dirname)
    file_iterator.iterate_over_files()

    mol_set_query = MoleculeSet.objects.filter(set_name=set_name)
    if mol_set_query:
        mol_set = mol_set_query[0]
    else:
        print("No such MoleculeSet!", file=sys.stderr)
        raise ValueError

    mol_set.molecules.add(*file_iterator.get_mol_list())
    mol_set.save()


class MoleculeIterator:
    """
    Base class for the iteration over all molecules in an SDF file
    """

    def __init__(self,
                 image_dir: str = os.path.join(settings.MEDIA_ROOT, 'images'),
                 save_model_list: bool = True):
        """
        :param image_dir: dir_path to the directory, in which the images
            will be saved
        :param save_model_list: if True, Molecule instances will be saved and
            can later be added to the set
        """

        self.save_model_list = save_model_list
        self.image_dir = image_dir
        self.molname_search_list = (
            'PUBCHEM_IUPAC_NAME', 'chembl_pref_name', 'GENERIC_NAME', "_Name")
        self.mol_list = []
        self.err_msgs = []

        if not os.path.isdir(self.image_dir):
            os.makedirs(self.image_dir)

    def iterate_over_molecules(self, file_path: str) -> bool:
        """
        Function for the iteration over all molecules in the provided filename;
        molecules will be saved in the Molecule model
        :param file_path: dir_path to an SDF file
        :return: True if the file is valid, else False
        """
        supplier = self.__get_supplier(file_path)
        if supplier is None:
            return False

        with supplier(file_path) as suppl:
            for mol in suppl:
                if not mol:
                    self._err_msg("Cannot parse molecule")
                    continue

                mol_props_instance = MoleculeProperties(mol, self.image_dir)
                mol_props_instance.search_for_name(self.molname_search_list)
                mol_instance = mol_props_instance.save_molecule()

                if self.save_model_list:
                    self.mol_list.append(mol_instance)
        return True

    def __get_supplier(self, file_path: str):
        """
        Validates file exists and has a valid extension and returns a matching
        Supplier
        :param file_path: dir_path to an SDF file
        :return: Supplier for parsing the file.
            If file is invalid, returns None
        """
        if not os.path.isfile(file_path):
            self._err_msg(f"{file_path} does not exist")
            return None

        if file_path.lower().endswith('.sdf'):
            return Chem.SDMolSupplier
        else:
            self._err_msg(f"{file_path} is not an SDF file! ")
            return None

    def get_mol_list(self) -> list[Molecule]:
        """
        Returns a list of all Molecule instances which have already
        been processed
        :return: list of Molecule instances
        """
        return self.mol_list

    def get_err_msgs(self) -> list[str]:
        """
        Returns a list of errors which occurred during the iteration
        :return: list of error messages
        """
        return self.err_msgs

    def _err_msg(self, err_msg: str) -> None:
        """
        Prints and logs the given error message
        :param err_msg: error that occurred
        """
        self.err_msgs.append(err_msg)
        print(err_msg, file=sys.stderr)

    def add_mol_list_to_set(self, molecule_set: MoleculeSet) -> None:
        """
        Insert all molecules into the given MoleculeSet
        """
        molecule_set.molecules.add(*self.mol_list)
        molecule_set.save()

    def add_to_set(self, set_name: str, set_description: Optional[str]) -> None:
        """
        Creates a new set with the given set_name and set_description if the
        set does not already exist and adds links the Molecules to the set
        :param set_name: name of the set
        :param set_description: description of the set
        """
        if set_description:
            mol_set = MoleculeSet(
                set_name=set_name,
                description=set_description)
        else:
            mol_set = MoleculeSet(
                set_name=set_name)

        mol_set.save()
        self.add_mol_list_to_set(mol_set)


class FileIterator(MoleculeIterator):
    """
    Class for the iteration over all files in a directory,
    adding the molecules into the DB.
    """

    def __init__(self, dirname: str):
        """
        :param dirname: dir_path of the directory
        """
        super().__init__()
        self.dirname = dirname

        if not os.path.isdir(dirname):
            err = f"{dirname} is not a valid directory"
            self._err_msg(err)
            raise ValueError(err)

    def iterate_over_files(self) -> None:
        """
        Function for the iteration over files in the given directory.
        Molecules will be saved in the Molecule model
        """
        for root, dirs, files in os.walk(self.dirname):
            for filename in files:
                path_to_file = os.path.join(root, filename)
                self.iterate_over_molecules(path_to_file)
