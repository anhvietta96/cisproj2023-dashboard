from compounds.calc_mol_properties import insert_into_existing_set
from dashboard.settings import MEDIA_ROOT


def run(set_id) -> None:
    """
    Adds the compounds from MEDIA_ROOT into the DB and links them to an
    already existing set with id set_id
    If set does not exist, raise ValueError
    """
    set_id = int(set_id)
    insert_into_existing_set(MEDIA_ROOT, set_id)
