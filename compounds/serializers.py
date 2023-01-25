from rest_framework import serializers
from .models import Molecule


class MoleculeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Molecule
        fields = ('inchi_key', 'log_p', 'molecular_formula', 'molecular_weight',
                  'num_h_acceptors', 'num_h_donors', 'num_rotatable_bonds',
                  'num_rings')
