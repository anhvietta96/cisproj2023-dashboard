import os.path
from django.db import models
from django.core.validators import MinLengthValidator, MinValueValidator
from dashboard.settings import MEDIA_ROOT


class testMolecule(models.Model):
    inchi_key = models.CharField(
        primary_key=True, max_length=27, validators=[MinLengthValidator(27)])

    log_p = models.FloatField(
        default=None, blank=True, null=True)
    molecular_formula = models.CharField(
        default=None, blank=True, null=True, max_length=50,
        validators=[MinLengthValidator(1)])
    molecular_weight = models.FloatField(
        default=None, blank=True, null=True, validators=[MinValueValidator(0)])
    num_h_acceptors = models.PositiveSmallIntegerField(
        default=None, blank=True, null=True)
    num_h_donors = models.PositiveSmallIntegerField(

        default=None, blank=True, null=True)
    num_rotatable_bonds = models.PositiveSmallIntegerField(
        default=None, blank=True, null=True)
    num_rings = models.PositiveSmallIntegerField(
        default=None, blank=True, null=True)
    image = models.ImageField(
        default=None, blank=True, null=True, max_length=255,
        upload_to='media/images')

    def num_satisfied_lipinski_rules(self):
        lipinskis_ro5 = (self.molecular_weight < 500,
                         self.num_h_acceptors <= 10,
                         self.num_h_donors <= 5,
                         self.log_p < 5)

        return sum(lipinskis_ro5)

    def __str__(self):
        return f"inchi_key={self.inchi_key}"


class testMoleculeSet(models.Model):
    set_id = models.AutoField(
        primary_key=True)
    set_name = models.CharField(
        max_length=120, validators=[MinLengthValidator(1)])
    description = models.CharField(
        max_length=250, default=None, blank=True, null=True)
    molecules = models.ManyToManyField(testMolecule)

    class Meta:
        ordering = ['set_name']

    def save(self, *args, **kwargs):
        super().save(*args, **kwargs)
        if not self.set_name:
            self.set_name = f"Set {self.set_id}"
            super().save(update_fields=['set_name'])

    def __str__(self):
        return self.set_name
