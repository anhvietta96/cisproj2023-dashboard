import dataclasses
from django.shortcuts import render
from django.views.generic import ListView
from .serializers import MoleculeSerializer
from .models import Molecule, MoleculeSet
from rest_framework import viewsets
from django.http import HttpRequest
from django.core.exceptions import ObjectDoesNotExist


@dataclasses.dataclass
class Set:
    set_name: str
    mol_list: list


def main_compound_view(request: HttpRequest):
    set_list = [Set(mol_set.set_name, mol_set.molecules.all())
                for mol_set in MoleculeSet.objects.all()]
    if request.method == 'POST':
        # Apply filter on molecules here
        pass
    return render(request, 'compounds/compounds.html', {'set_list': set_list})


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
    q = Molecule.objects.filter(name__icontains=query)
    search_results = []
    property_list = Molecule.objects.get_all_attr()

    property_list.remove('image')
    property_list = ['image'] + property_list
    display_property_list = [property.replace('_', ' ').title() for property in
                             property_list]

    for mol in q:
        value_list = []
        for property in property_list:
            if property != 'image':
                value_list.append(getattr(mol, property))
            else:
                value_list.append(mol.image.url)
        search_results.append(value_list)

    data = {'header': display_property_list, 'table': search_results}
    return render(request, 'search.html', data)
