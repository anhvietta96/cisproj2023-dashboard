from rest_framework import serializers
from .models import Substance


class SubstanceSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Substance
        fields = ('PubChem_Name', 'PubChem_Alias', 'PubChem_CID', 'PubChem_Mass', 'PubChem_Formula')
