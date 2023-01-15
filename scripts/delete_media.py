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
        print(f"Folder {dir_path} does not exist!", file=sys.stderr)
    os.mkdir(dir_path)


def run(*args: str) -> None:
    """
    Deletes (recursively) the content of MEDIA_ROOT
    Deletes the content of MEDIA_ROOT/uploaded_data/
    :param args: list of arguments:
                 sdf: delete the content of MEDIA_ROOT/uploaded_data/ instead
    """
    if "sdf" in args:
        path = os.path.join(MEDIA_ROOT, 'uploaded_data/')
        delete_dir_content(path)
    else:
        delete_dir_content(MEDIA_ROOT)
