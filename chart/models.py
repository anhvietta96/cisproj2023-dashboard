from django.db import models
from compounds.models import Molecule
# Create your models here.

properties = Molecule.__dict__["__doc__"]
property_list = properties[9:].replace(',','').replace(')','').split()

PROPERTY_CHOICES = tuple((property,property) for property in property_list)

class PropertyChoice(models.Model):
    property = models.CharField(max_length=20,choices=PROPERTY_CHOICES)