"""
Methods to calculate molecular properties and draw molecule images
"""
import os.path
from .models import Molecule
from typing import Optional, Iterable
from django.conf import settings
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, inchi, \
    Crippen, Draw, rdchem


class MoleculeProperties:
    """
    Class to calculate molecular properties and draw molecule images
    """

    def __init__(self, mol: rdchem.Mol, rel_image_dir: Optional[str] = None):
        """
        :param mol: molecule of interest, provided by a supplier
        :param rel_image_dir: relative dir_path (from MEDIA_ROOT) to directory
            in which the generated image will be saved
        """
        self.mol = mol
        self.rel_image_dir = rel_image_dir
        self.name = None

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

    def search_for_name(self, name_list: Iterable[str]) -> Optional[str]:
        """
        Searches for the name of the molecule in name_list and initializes
        self.name
        :param name_list: Iterable which includes list of properties describing
            the name of the molecule which could occur in self.mol
        :return: name of the molecule if the name has been found
        """
        for name_prop in name_list:
            if self.mol.HasProp(name_prop):
                name = self.mol.GetProp(name_prop)

                if name and name != 'None':
                    self.name = name
                    return self.name

    def draw_image(self) -> Optional[str]:
        """
        Draws an image of the molecule and saves it to the given
        image directory.
        :return: dir_path to the image file if self.rel_image_dir is not None
        """
        if self.rel_image_dir is not None:
            image_path = os.path.join(settings.MEDIA_ROOT,
                                      self.rel_image_dir,
                                      f"{self.inchi_key}.png")
            Draw.MolToFile(self.mol, image_path)
            return image_path

    def return_dict(self) -> dict:
        """
        Returns a dictionary of molecule properties
        :return: dictionary of molecular descriptors
        """
        return {
            'inchi_key': self.inchi_key,
            'name': self.name,
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
        Saves molecular properties into the DB
        :return: Instance of saved Molecule
        """
        image_path = self.draw_image()

        m = Molecule(
            inchi_key=self.inchi_key,
            name=self.name,
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
