from rest_framework import serializers
from .models import Molecule

class MoleculeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Molecule
        fields = ('inchi_key','molecular_mass','num_h_acceptors','num_h_donors')