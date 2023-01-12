from django.db import models
from django.core.validators import MinLengthValidator, MinValueValidator

class MoleculeManager(models.Manager):
    def get_all_attr(self):
        all_attr_list = []
        for field in self.model._meta.get_fields():
            if 'Field' in str(type(self.model._meta.get_field(field.name))):
                all_attr_list.append(field.name)
        return all_attr_list

    def get_num_attr(self):
        num_attr_list = ['Integer','Float']
        num_field = []
        all_field_name = [field.name for field in self.model._meta.get_fields()]
        for field_name in all_field_name:
            for attr in num_attr_list:
                if attr in str(type(self.model._meta.get_field(field_name))):
                    num_field.append(field_name)
        return num_field

    def get_str_attr(self):
        num_attr = self.get_num_attr()
        all_field_name = [field.name for field in self.model._meta.get_fields()]
        return [attr for attr in all_field_name and attr not in num_attr]


class Molecule(models.Model):
    objects = MoleculeManager()


    inchi_key = models.CharField(
        primary_key=True, max_length=27, validators=[MinLengthValidator(27)])
    name = models.TextField(
        default=None, blank=True, null=True)

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
        upload_to='images/')

    def num_satisfied_Lipinski_rules(self) -> int:
        """
        Returns how many rules of Lipinski's rule of five are fulfilled
        :return: number of satisfied rules
        """
        lipinskis_ro5 = (self.molecular_weight < 500,
                         self.num_h_acceptors <= 10,
                         self.num_h_donors <= 5,
                         self.log_p < 5)

        return sum(lipinskis_ro5)

    def __str__(self):
        return self.inchi_key


class MoleculeSet(models.Model):
    set_name = models.CharField(
        primary_key=True, max_length=120, validators=[MinLengthValidator(1)])
    description = models.CharField(
        max_length=250, default=None, blank=True, null=True)
    molecules = models.ManyToManyField(Molecule)

    class Meta:
        ordering = ['set_name']

    def __str__(self):
        return self.set_name
