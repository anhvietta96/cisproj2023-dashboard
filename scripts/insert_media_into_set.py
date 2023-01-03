from typing import Optional
from compounds.file_handler import FileIterator
from compounds.models import MoleculeSet
from dashboard.settings import MEDIA_ROOT


def insert_into_set(dirname: str, set_name: str) -> None:
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


def run(set_name: str, dir_name: Optional[str] = None) -> None:
    """
    Adds the compounds from MEDIA_ROOT into the DB and links them to an
    set with name set_name
    :param set_name: name of the set
    :param dir_name: dir_path to the directory of SDF files.
        If dir_name is None, the directory is  MEDIA_ROOT
    """
    if dir_name is None:
        insert_into_set(MEDIA_ROOT, set_name)
    else:
        insert_into_set(dir_name, set_name)

