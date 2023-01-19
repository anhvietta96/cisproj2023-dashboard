from django.test import SimpleTestCase
from django.urls import reverse, resolve
from compounds.views import main_compound_view, molecule_single_view, SearchResultsView, CompoundViewSet


class TestUrls(SimpleTestCase):

    def test_compounds_url_resolves(self):
        url = reverse('compounds')
        #print(resolve(url))
        self.assertEquals(resolve(url).func, main_compound_view)

    #def test_search_results_view(self):
     #   url = reverse('search_url')
      #  #print(resolve(url))
       # self.assertEquals(resolve(url).func, SearchResultsView) #here one error and one failure

    def test_search_results_view(self):
        url = '/search/'
        self.assertEquals(resolve(url).view_name, 'search_url') #here two errors
