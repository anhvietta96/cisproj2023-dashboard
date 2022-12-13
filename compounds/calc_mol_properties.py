"""
Methods for the calculation of molecular properties
"""
import os.path
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, inchi, \
    Crippen, Draw
from .models import Molecule


class MoleculeProperties:
    """
    Class to calculate molecular properties
    """

    def __init__(self, mol, image_dir: str):
        """
        :param mol: molecule of interest, provided by a supplier
        """
        self.mol = mol
        self.image_dir = image_dir

        # primary key
        self.inchi_key = inchi.MolToInchiKey(self.mol)

        # properties for Lipinksi's rule of five
        self.log_p = round(Crippen.MolLogP(self.mol), 6)
        self.molecular_weight = round(Descriptors.MolWt(self.mol), 4)
        self.num_h_acceptors = Lipinski.NumHAcceptors(self.mol)
        self.num_h_donors = Lipinski.NumHDonors(self.mol)

        # additional properties
        self.molecular_formula = rdMolDescriptors.CalcMolFormula(self.mol)
        self.num_rings = rdMolDescriptors.CalcNumRings(self.mol)
        self.rotatable_bonds = Lipinski.NumRotatableBonds(self.mol)

    def draw_image(self) -> str:
        """
        Draws an image of the molecule and saves it to an existing directory.
        :return: path to the image file
        """
        image_path = os.path.join(self.image_dir, f"{self.inchi_key}.png")
        Draw.MolToFile(self.mol, image_path)
        return image_path

    def return_dict(self) -> dict:
        """
        Returns a dictionary of molecule properties
        :return: dictionary of molecular descriptors
        :rtype: dict
        """
        return {
            'inchi_key': self.inchi_key,
            'log_p': self.log_p,
            'molecular_formula': self.molecular_formula,
            'molecular_weight': self.molecular_weight,
            'num_h_acceptors': self.num_h_acceptors,
            'num_h_donors': self.num_h_donors,
            'num_rings': self.num_rings,
            'num_rotatable_bonds': self.rotatable_bonds
        }

    def save_molecule(self) -> Molecule:
        """
        Saves molecular properties
        :return: Instance of saved Molecule
        :rtype: Molecule
        """
        image_path = self.draw_image() if self.image_dir else None

        m = Molecule(
            inchi_key=self.inchi_key,
            log_p=self.log_p,
            molecular_formula=self.molecular_formula,
            molecular_weight=self.molecular_weight,
            num_h_acceptors=self.num_h_acceptors,
            num_h_donors=self.num_h_donors,
            num_rings=self.num_rings,
            num_rotatable_bonds=self.rotatable_bonds,
            image=image_path
        )
        m.save()
        return m

    def __str__(self):
        return str(self.return_dict())
