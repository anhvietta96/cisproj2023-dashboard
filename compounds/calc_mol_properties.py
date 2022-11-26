from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, inchi
from .models import Molecule


class MoleculeIteratorClass:
    def __init__(self, filename: str, supplier=Chem.SDMolSupplier):
        self.filename = filename
        self.supplier = supplier
        self.mol_list = []

    def iterate_over_molecules(self):
        with self.supplier(self.filename) as suppl:
            for mol in suppl:
                if not mol:
                    continue

                mol_props = MoleculeProperties(mol)
                mol_props.save_molecule()
                self.mol_list.append(mol_props)

    def printall(self):
        for mol in self.mol_list:
            print(mol)


class MoleculeProperties:
    def __init__(self, mol):
        self.mol = mol

        # primary key
        self.inchi_key = inchi.MolToInchiKey(self.mol)

        # properties for Lipinksi's rule of five
        self.molecular_weight = round(Descriptors.MolWt(self.mol), 4)
        self.log_p = round(Descriptors.MolLogP(self.mol), 6)
        self.num_h_acceptors = Lipinski.NumHAcceptors(self.mol)
        self.num_h_donors = Lipinski.NumHDonors(self.mol)

        # additional properties
        self.rotatable_bonds = Lipinski.NumRotatableBonds(self.mol)
        self.molecule_formula = rdMolDescriptors.CalcMolFormula(self.mol)

    def return_dict(self):
        return {
            'inchi_key': self.inchi_key,
            'molecular_weight': self.molecular_weight,
            'log_p': self.log_p,
            'rotatable_bonds': self.rotatable_bonds,
            'num_h_acceptors': self.num_h_acceptors,
            'num_h_donors': self.num_h_donors,
            'molecule_formula': self.molecule_formula,
        }

    def save_molecule(self):
        m = Molecule(
            inchi_key=self.inchi_key,
            num_h_acceptors=self.num_h_acceptors,
            num_h_donors=self.num_h_donors,
            log_p=self.log_p,
            molecular_mass=self.molecular_weight
        )
        m.save()

    def __str__(self):
        return str(self.return_dict())


if __name__ == '__main__':
    filename = 'Beispieldateien/Drugs/Drugs.sdf'
    x = MoleculeIteratorClass(filename)
    x.iterate_over_molecules()
    x.printall()
