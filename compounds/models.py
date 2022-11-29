from django.db import models


class Molecule(models.Model):
    inchi_key = models.CharField(primary_key=True, max_length=27)

    log_p = models.FloatField(
        default=None, blank=True, null=True)
    num_h_acceptors = models.PositiveSmallIntegerField(
        default=None, blank=True, null=True)
    num_h_donors = models.PositiveSmallIntegerField(
        default=None, blank=True, null=True)
    molecular_mass = models.FloatField(
        default=None, blank=True, null=True)

    def __str__(self):
        return f"inchi_key={self.inchi_key}"


class MoleculeSet(models.Model):
    set_id = models.AutoField(primary_key=True)
    set_name = models.CharField(max_length=120)
    molecules = models.ManyToManyField(Molecule)

    class Meta:
        ordering = ['set_name']

    def __str__(self):
        return self.set_name
