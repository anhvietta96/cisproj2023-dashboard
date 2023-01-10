from django.shortcuts import render
from django.views.generic import ListView
from .serializers import MoleculeSerializer
from .models import Molecule, MoleculeSet
from rest_framework import viewsets
from django.http import HttpRequest
from django.core.exceptions import ObjectDoesNotExist


def main_compound_view(request: HttpRequest):
    sets = MoleculeSet.objects.all()
    return render(request, 'compounds/compounds.html', {'sets': sets})


class CompoundViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all().order_by('inchi_key')
    serializer_class = MoleculeSerializer
    http_method_names = ['get', 'head']


def molecule_single_view(request: HttpRequest, inchi_key: str):
    try:
        molecule = Molecule.objects.get(pk=inchi_key)
    except ObjectDoesNotExist:
        return render(request, 'compounds/compounds.html')

    return render(request, 'compounds/compounds.html',
                  {'molecules': [molecule]})


class SearchResultsView(ListView):
    model = Molecule
    template_name = 'search.html'

    def get_queryset(self):
        query = self.request.GET.get("search_query")
        search_result = Molecule.objects.filter(inchi_key__icontains=query)
        return search_result


def search_results(request):
    query = request.GET.get("search_query")
    q = Molecule.objects.filter(inchi_key__icontains=query)
    search_results = []
    properties = Molecule.__dict__["__doc__"]
    property_list = properties[9:].replace(',', '').replace(')', '').split()
    
    for mol in q:
        value_list = []
        for property in property_list:
            if property != 'image':
                value_list.append(getattr(mol, property))
            else:
                value_list.append(mol.image.url)
        search_results.append(value_list)
    
    property_list.remove('image')
    property_list = ['image'] + property_list

    data = {'header': property_list, 'table': search_results}
    return render(request, 'search.html', data)
