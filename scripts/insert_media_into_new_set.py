from compounds.calc_mol_properties import insert_into_new_set
from dashboard.settings import MEDIA_ROOT


def run(set_name: str = '') -> None:
    """
    Adds the compounds from MEDIA_ROOT into the DB and links them to a new set
    with name set_name
    """
    insert_into_new_set(MEDIA_ROOT, set_name)
