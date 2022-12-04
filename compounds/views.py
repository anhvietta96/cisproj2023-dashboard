from django.shortcuts import render
import os
from glob import glob
from pathlib import Path
from .moleculeClass import moleculeClass
from dashboard.settings import MEDIA_ROOT
from django.views.generic import ListView
from .serializers import MoleculeSerializer
from .models import Molecule
from rest_framework import viewsets

def main_compound_view(request):
    properties = ['PUBCHEM_COMPOUND_CID',
                  'PUBCHEM_CACTVS_HBOND_ACCEPTOR',
                  'PUBCHEM_CACTVS_HBOND_DONOR',
                  'PUBCHEM_IUPAC_INCHI',
                  'PUBCHEM_MOLECULAR_WEIGHT',
                  'PUBCHEM_MOLECULAR_FORMULA',
                  'PUBCHEM_XLOGP3']

    mol_properties = []
    mol_names = [Path(file).stem for file in
                 glob(os.path.join(MEDIA_ROOT, '*.sdf'))]

    for mol_name in mol_names:
        mol_instance = moleculeClass("media/", mol_name)
        mol_instance.loadMolecule()
        for property_dict in mol_instance.getPropDictIterator():
            temp = [mol_name]
            for prop in properties:
                temp.append((prop, property_dict.get(prop)))
            mol_properties.append(temp)

    return render(request, 'compounds/compounds.html',
                  {'mols_with_properties': mol_properties})

class CompoundViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all().order_by('inchi_key')
    serializer_class = MoleculeSerializer
    http_method_names = ['get','head']

class SearchResultsView(ListView):
    model = Molecule
    template_name = 'search.html'

    def get_queryset(self):
        query=self.request.GET.get("search_query")
        search_result = Molecule.objects.filter(inchi_key__icontains=query)
        return search_result

def search_results(request):
    query=request.GET.get("search_query")
    q = Molecule.objects.filter(inchi_key__icontains=query)
    search_results = []
    properties = Molecule.__dict__["__doc__"]
    property_list = properties[9:].replace(',','').replace(')','').split()
    for mol in q:
        value_list = []
        for property in property_list:
            value_list.append(getattr(mol,property))
        search_results.append(value_list)
    data = {'table':search_results}
    return render(request,'search.html',data)