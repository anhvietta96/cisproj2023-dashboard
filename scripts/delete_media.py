import os.path
import shutil
import sys
from dashboard.settings import MEDIA_ROOT


def delete_dir_content(path):
    try:
        shutil.rmtree(path)
    except FileNotFoundError:
        print("Folder does not exists!", file=sys.stderr)
    os.mkdir(path)


def run(*args) -> None:
    """
    Delete the content of MEDIA_ROOT/uploaded_data/(recursively)
    :param args: list of arguments:
                 all: delete the content of MEDIA_ROOT instead
    """
    if 'all' in args:
        delete_dir_content(MEDIA_ROOT)

    else:
        path = os.path.join(MEDIA_ROOT, 'uploaded_data/')
        delete_dir_content(path)
