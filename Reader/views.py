from rest_framework import viewsets
from .serializers import SubstanceSerializer
from .models import Substance


class SubstanceViewSet(viewsets.ModelViewSet):
    queryset = Substance.objects.all().order_by('PubChem_CID')
    serializer_class = SubstanceSerializer
    http_method_names = ['get', 'head']
