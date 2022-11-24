from __future__ import print_function
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import Draw
import glob
from rdkit.Chem import PandasTools
import pandas as pd
import os
from rdkit import RDConfig

class moleculeClass():
    """
    First part of implementation of RDKit-features
    """
    def __init__(self, dir, name):
        self.dir = dir
        self.name = name

    def buildLink(self):
        """
        This function builds the directory link to
        the related sdf file.
        :param: dir: Link to directory , name: name of the file
        :return: Link to the .sdf-File
        :rtype: string
        """
        #assert self.dir and self.name is not none
        return self.dir+self.name+".sdf"

    def loadMolecule(self):
        """
        Loads sdf-Files.
        :param: dir: Linkt to directory, name: name of the file
        :return: Supplier 
        :rtype: Supplier
        """
        self.link = self.buildLink()
        self.suppl= Chem.SDMolSupplier(self.link)

        #Wird nicht mehr benötigt!
    def getHeavyAtomNum(self):
        """
        Number of heavy atoms.
        :param: suppl: Supplier
        :return: Number of heavy atoms.
        :rtype: float 
        """
        for mol in self.suppl:
            if mol is None: continue #test molecule before using it#
            list=mol.GetNumAtoms()
        return list

    def drawMolecule2D(self):
        """
        Provides a picture of the actual supplier
        :param: suppl: Supplier, name: name of file to be stored
        :return: none
        :rtype: none
        """
        ms = [x for x in self.suppl if x is not None]
        for m in ms: tmp=AllChem.Compute2DCoords(m)
        Draw.MolToFile(ms[0], self.name+".png")

    def molToSmiles(self):
        """
        Mol to smiles
        :param: mol: Molecule to be converted
        :return: Smiles-Type
        :rtype: Smiles
        """
        smi = Chem.MolToSmiles(self.mol)
        return smi

    def getSubStructureNum(self,struct):
        """
        Get num of matching substructure
        :param: suppl: Supplier, struct: pattern
        :return: number of matches
        :rtype: Integer
        """
        patt=Chem.MolFromSmarts(struct)
        matches = [m for m in self.suppl if m.HasSubstructMatch(patt)]
        return len(matches)
    #_______________________________________________________________________________#

    #Ab hier beginnen die Eigenschaften#


    def getPropDictIterator(self):
        """
        Get a dictonary of atomar discriptors
        :param: mol: Molecul of interest
        :return: dictonary with all possible discriptors
        :rtype: Dictionary
        """
        for mol in self.suppl:
            yield mol.GetPropsAsDict()


    #Lipinski Rule of 5
    """
    Descriptors.MolLogP(mol)

    Descriptors.ExactMolWt(mol)

    Descriptors.NumHDonors(mol)
    Descriptors.NumHAcceptors(mol)
    #Others
    Descriptors.NumRotatableBonds(mol)
    Descriptors.RingCount(mol)
    Descriptors.BalabanJ(mol) #Erlaubt detailreiche Unterscheidung ähnlicher Strukturen
    Descriptors.BertzCT(mol) #Komplexitaet vergleichen
    """

        #Pandas
"""
frame = PandasTools.LoadSDF('Arginine.sdf',smilesName='SMILES',molColName='Molecule',
           includeFingerprints=True)
frame.values
frame.columns
"""