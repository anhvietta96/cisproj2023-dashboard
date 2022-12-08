from django.shortcuts import render
from django.views.generic.edit import CreateView
from compounds.models import Molecule
'''from .models import PropertyChoice
from .forms import PropertyForm'''
import json

class Chart:
    pass

# Create your views here.
def ChartOptions(request):
    data = {}
    properties = Molecule.__dict__["__doc__"]
    property_list = properties[9:].replace(',','').replace(')','').split()
    data['property_list']=property_list
    return render(request,'chart_options.html',data)

'''
def ChartResult(request):
    
    x_axis = int(request.POST.get('x-axis')[0])-1
    y_axis = int(request.POST.get('y-axis')[0])-1
    q = Molecule.objects.all()
    
    properties = Molecule.__dict__["__doc__"]
    property_list = properties[9:].replace(',','').replace(')','').split()

    ChartOptions = {"chart":{"height":350,"type":"scatter","zoom":{"enabled":True,"type":"xy"}},
    "xaxis": {
    "tickAmount": 10,"labels" : {"show":True}, "title":{"text":property_list[x_axis]}
    },"yaxis": {
    "tickAmount": 7,"labels" : {"show":True}, "title":{"text":property_list[y_axis]}}}

    data = []
    info = {"name":property_list[x_axis]+"/"+property_list[y_axis]}
    info["data"] = []
    for mol in q:
        info["data"].append([getattr(mol,property_list[x_axis]),getattr(mol,property_list[y_axis])])
    data.append(info)
    
    ChartOptions["series"]=data
    result = {'title':'Created Chart','options':ChartOptions}
    return render(request,'chart_result.html',result)
'''

def ChartResult(request):
    
    x_axis = int(request.POST.get('x-axis')[0])-1
    y_axis = int(request.POST.get('y-axis')[0])-1
    q = Molecule.objects.all()
    
    properties = Molecule.__dict__["__doc__"]
    property_list = properties[9:].replace(',','').replace(')','').split()

    ChartOptions = {}

    data = []
    for mol in q:
        mol_inf = []
        mol_inf.append(getattr(mol,property_list[x_axis]))
        mol_inf.append(getattr(mol,property_list[y_axis]))
        for property in property_list:
            if property != property_list[x_axis] and property != property_list[y_axis]:
                mol_inf.append(getattr(mol,property))
        data.append(mol_inf)

    ChartOptions['legend']=['Group 1']
    ChartOptions['header']=[property_list[x_axis],property_list[y_axis]]
    for property in property_list:
        if property != property_list[x_axis] and property != property_list[y_axis]:
            ChartOptions['header'].append(property)
    ChartOptions['name']=['Group 1']
    ChartOptions["data"]=[data]
    ChartOptions['title']=property_list[x_axis]+"/"+property_list[y_axis]
    result = {'title':'Created Chart','options':ChartOptions}
    
    return render(request,'chart_result.html',result)