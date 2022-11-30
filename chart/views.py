from django.shortcuts import render
from compounds.models import Molecule
import json

class Chart:
    pass

# Create your views here.
def ChartOptions(request):
    pass

def ChartResult(request):
    attr1 = 2
    attr2 = 3
    
    q = Molecule.objects.all()
    
    attributes = Molecule.__dict__["__doc__"]
    attributes_list = attributes[9:].replace(',','').replace(')','').split()

    ChartOptions = {"chart":{"height":350,"type":"scatter","zoom":{"enabled":True,"type":"xy"}},
    "xaxis": {
    "tickAmount": 10,"label" : attributes_list[attr1]
    },"yaxis": {
    "tickAmount": 7,"label" : attributes_list[attr2]}}

    data = []
    info = {"name":attributes_list[attr1]+"/"+attributes_list[attr2]}
    info["data"] = []
    for mol in q:
        info["data"].append([getattr(mol,attributes_list[attr1]),getattr(mol,attributes_list[attr2])])
    data.append(info)
    
    ChartOptions["series"]=data
    result = {'title':'Created Chart','options':json.dumps(ChartOptions)}
    print(result)
    return render(request,'chart.html',result)