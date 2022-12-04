from compounds.calc_mol_properties import FileIterator
from dashboard.settings import MEDIA_ROOT


def run(*args):
    """
    Adds the compounds from the directories in args into the DB
    Default directory is MEDIA_ROOT
    """
    if not args:
        args = [MEDIA_ROOT]

    for dirname in args:
        file_iterator = FileIterator(dirname)
        file_iterator.iterate_over_files()
