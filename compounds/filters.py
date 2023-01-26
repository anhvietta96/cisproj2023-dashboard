"""
        This function filters a Molecule-Query-Set
        :param:
                filter_query: one of    -->tuples [l,u] for numerical attributes, l being the lower and u being the upper bound.
                                        -->a String for attributes stored in charFields
                attribute: Name of the ModelField holding the attribute one wants to filter by


        :return: result
        :rtype: QuerySet
        """
from compounds.models import MoleculeSet


def filter_sets(attribute, filter_query):
    sets = MoleculeSet.objects.all()
    result = []
    for set in sets:
        result.append((set.molecules.filter(**build_molecules_filter(attribute, filter_query)), set.set_name))
    return result


def build_molecules_filter(attribute, filter_query):
    if type(filter_query) == list or type(filter_query) == tuple:
        return {
            f"{attribute}__gte": min(filter_query),
            f"{attribute}__lte": max(filter_query)
        }
    elif type(filter_query) == str:
        return {f"{attribute}__icontains": filter_query}
    else:
        raise TypeError(
            "Function takes either List of lower and upper bounds for numerical attributes or string for string attributes")
