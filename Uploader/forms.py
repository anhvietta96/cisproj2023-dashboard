from django import forms
from django.core.validators import FileExtensionValidator
from .models import SDFile


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
    set = forms.CharField(max_length=120)
    description = forms.CharField(max_length=250)
    documents = forms.FileField(
        widget=forms.ClearableFileInput(attrs={'multiple': True}),
        validators=[FileExtensionValidator(allowed_extensions=['sdf', 'smi'])])
