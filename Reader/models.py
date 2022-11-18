from django.db import models


class Substance(models.Model):
    PubChem_Name = models.CharField(max_length=50)
    PubChem_Alias = models.CharField(max_length=50)
    PubChem_CID = models.IntegerField(primary_key=True)
    PubChem_Mass = models.FloatField()
    PubChem_Formula = models.CharField(max_length=50)
    Local_Subcategory = models.TextChoices('Subgroup', 'DRUGS PROTEINOGENICAMINOACIDS')
