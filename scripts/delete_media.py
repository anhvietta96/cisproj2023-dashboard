import os.path
import shutil
import sys
from dashboard.settings import MEDIA_ROOT


def delete_dir_content(dir_path: str) -> None:
    """
    Deletes the given directory recursively
    :param dir_path: path to the directory
    """
    try:
        shutil.rmtree(dir_path)
    except FileNotFoundError:
        print(f"Folder {dir_path} does not exists!", file=sys.stderr)
    os.mkdir(dir_path)


def run(*args: str) -> None:
    """
    Deletes (recursively) the content of MEDIA_ROOT/uploaded_data/
    :param args: list of arguments:
                 all: delete the content of MEDIA_ROOT instead
    """
    if 'all' in args:
        delete_dir_content(MEDIA_ROOT)
    else:
        path = os.path.join(MEDIA_ROOT, 'uploaded_data/')
        delete_dir_content(path)
