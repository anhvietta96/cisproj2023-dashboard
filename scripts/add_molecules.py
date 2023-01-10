import os.path
import sys
from compounds.file_handler import FileIterator
from dashboard.settings import MEDIA_ROOT


def run(*args: str) -> None:
    """
    Adds the compounds from the directories in args into the DB
    Default directory is MEDIA_ROOT
    :param args: list of directory paths
    """
    if not args:
        args = [MEDIA_ROOT]

    for dirname in args:
        if not os.path.isdir(dirname):
            print(f"{dirname} is not a valid directory", file=sys.stderr)
            continue

        file_iterator = FileIterator(dirname)
        file_iterator.iterate_over_files()
