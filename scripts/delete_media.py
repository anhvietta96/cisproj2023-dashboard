import shutil
from dashboard.settings import MEDIA_ROOT


def run() -> None:
    """
    Deletes MEDIA_ROOT (recursively)
    """
    shutil.rmtree(MEDIA_ROOT)
