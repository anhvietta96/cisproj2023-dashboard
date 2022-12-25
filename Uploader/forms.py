from django import forms
from django.core.validators import FileExtensionValidator
from .models import SDFile
from compounds.models import MoleculeSet


# Single file upload only, based on a model so name conflicts get automatically
# resolved
class SDFileForm(forms.ModelForm):
    class Meta:
        model = SDFile
        fields = ('description', 'document')
        widget = \
            {'document': forms.ClearableFileInput(attrs={'multiple': True})}


# Multiple file upload, must write resolution for name conflicts
class SDFileMult(forms.Form):
    set_name = forms.CharField(max_length=120)
    set_already_exists = forms.BooleanField(required=False)
    description = forms.CharField(max_length=250, required=False)
    documents = forms.FileField(
        widget=forms.ClearableFileInput(attrs={'multiple': True}),
        validators=[FileExtensionValidator(allowed_extensions=['sdf', 'smi'])])

    def clean(self):
        cleaned_data = super().clean()
        already_exists = cleaned_data["set_already_exists"]
        set_name = cleaned_data["set_name"]

        mol_set_query = MoleculeSet.objects.filter(set_name=set_name)
        if mol_set_query:
            if already_exists is False:
                raise forms.ValidationError(
                    f"A set with the name {set_name} already exists!")
        else:
            if already_exists is True:
                raise forms.ValidationError(
                    f"A set with the name {set_name} does not exist!")

        return cleaned_data

