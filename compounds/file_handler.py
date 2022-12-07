"""
File handling classes and functions for .sdf (and .smi) files
"""
import os.path
import sys
from rdkit import Chem
from .models import Molecule, MoleculeSet
from .calc_mol_properties import MoleculeProperties
from dashboard.settings import MEDIA_ROOT


def insert_into_new_set(dirname: str, set_name: str = '') -> None:
    """
    Adds the compounds from dirname into the DB and links them to a new
    MoleculeSet with name set_name
    :param dirname: name of the directory which includes SDF or SMI files
    :param set_name: name of the new MoleculeSet
    """
    file_iterator = FileIterator(dirname)
    file_iterator.iterate_over_files()

    mol_set = MoleculeSet(set_name=set_name)
    mol_set.save()

    mol_set.molecules.add(*file_iterator.get_mol_list())
    mol_set.save()


def insert_into_existing_set(dirname: str, set_id: int = None) -> None:
    """
    Adds the compounds from dirname into the DB and links them to an already
    existing MoleculeSet with id set_id
    :param dirname: name of the directory which includes SDF or SMI files
    :param set_id: set_id of a specific MoleculeSet
    """
    file_iterator = FileIterator(dirname)
    file_iterator.iterate_over_files()

    mol_set_query = MoleculeSet.objects.filter(set_id=set_id)
    if mol_set_query:
        mol_set = mol_set_query[0]
    else:
        print("No such MoleculeSet!", file=sys.stderr)
        raise ValueError

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
        self.err_msgs = []

        if not os.path.isdir(dirname):
            print(f"{dirname} is not a valid directory", file=sys.stderr)
            raise ValueError

    def iterate_over_files(self) -> None:
        for root, dirs, files in os.walk(self.dirname):
            for filename in files:
                path_to_file = os.path.join(root, filename)
                self._handle_mol_iterator(path_to_file)

    def _handle_mol_iterator(self, path_to_file) -> bool:
        """
        :param path_to_file: path to a .sdf or .smi file
        :return: True if MoleculeIterator was successfully, else False
        """
        try:
            mol_iterator = MoleculeIterator(path_to_file)
        except ValueError:
            tmp_err_msg = f"Skipped {path_to_file}"
            self.err_msgs.append(tmp_err_msg)
            print(tmp_err_msg, file=sys.stderr)
            return False

        mol_iterator.iterate_over_molecules()
        self.mol_list.extend(mol_iterator.get_mol_list())
        return True

    def get_mol_list(self) -> list[Molecule]:
        """
        Returns a list of all Molecule instances from self.dirname
        """
        return self.mol_list

    def add_mol_list_to_set(self, molecule_set: MoleculeSet) -> None:
        """
        Insert all molecules into the given MoleculeSet
        """
        molecule_set.molecules.add(*self.mol_list)
        molecule_set.save()

    def get_err_msgs(self) -> list[str]:
        """
        Returns a list of error messages
        """
        return self.err_msgs


class MoleculeIterator:
    """
    Class for the iteration over all molecules in an .sdf or .smi file
    """

    def __init__(self,
                 filename: str,
                 image_dir: str = os.path.join(MEDIA_ROOT, 'images'),
                 save_model_instances: bool = True,
                 print_mol_dicts: bool = False,
                 multithreaded: bool = False):

        self.filename = filename
        self.supplier = None
        self.image_dir = image_dir
        self.multithreaded = multithreaded
        self.print_mol_dicts = print_mol_dicts
        self.save_model_list = save_model_instances
        self.mol_list = []
        self.err_msgs = []

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
        Function for the iteration over all molecules in the provided filename;
        molecules will be saved in the Molecule model
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
                    print(mol_props_instance)

                if self.image_dir is not None:
                    mol_props_instance.draw_image(self.image_dir)

    def get_mol_list(self) -> list[Molecule]:
        """
        Returns a list of all Molecule instances from self.filename
        """
        return self.mol_list
