from compounds.models import Molecule, MoleculeSet


def run():
    molecules = Molecule.objects.all()
    molecules.delete()

    molecule_sets = MoleculeSet.objects.all()
    molecule_sets.delete()
