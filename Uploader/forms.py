from django import forms
from .models import SDFile

#Single file upload only, based on model so name conficts get automatically resolved
class SDFileForm(forms.ModelForm):
    class Meta:
        model=SDFile
        fields=('description','document')
        widget={'document': forms.ClearableFileInput(attrs={'multiple':True})}

#Multiple file upload, must write resolution for name conflicts
class SDFileMult(forms.Form):
    set = forms.CharField(max_length=100)
    description = forms.CharField(max_length=100)
    documents = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple':True}))