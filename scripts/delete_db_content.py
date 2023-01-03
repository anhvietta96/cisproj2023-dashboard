from compounds.models import Molecule, MoleculeSet


def run(*args: str) -> None:
    """
    Deletes all Molecule and MoleculeSet instances in the DB
    :param args: list of arguments
                 no-sets: exclude all MoleculeSets
                 no-molecules: exclude all Molecules
    :return: None
    """
    if 'no-molecules' not in args:
        molecules = Molecule.objects.all()
        molecules.delete()

    if 'no-sets' not in args:
        molecule_sets = MoleculeSet.objects.all()
        molecule_sets.delete()
