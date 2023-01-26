from django.test import SimpleTestCase
from django.urls import reverse, resolve
from compounds.views import main_compound_view, molecule_single_view, SearchResultsView, CompoundViewSet


class TestUrls(SimpleTestCase):

    def test_search_results_view(self):
        url = '/search/'
        self.assertEquals(resolve(url).view_name, 'search_url')

    def test_1(self):
        self.assertEquals(1, 1)

    def test_2(self):
        self.assertEquals("Ringsssszzz", "Ringsssszzz")
