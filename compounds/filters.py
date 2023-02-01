"""
            This function filters a all MoleculeSets by the specified query by the specified attribute
            :param:
                    filter_query: one of    -->tuples [l,u] for numerical attributes, l being the lower and u being the upper bound.
                                            -->a String for attributes stored in charFields
                    attribute: Name of the ModelField holding the attribute one wants to filter by
                    set_listim:


            :return: result
            :rtype: list of tuples, each tuple contains the filtered set and the name of the set
            """

attr_name_mapping = {
    "Formula": "molecular_formula",
    "logP": "log_p",
    "Weight": "molecular_weight",
    "H-acceptors": "num_h_acceptors",
    "H-donors": "num_h_donors",
    "Rotatable bonds": "num_rotatable_bonds",
    "Rings": "num_rings"
}


def filter_sets(attribute, filter_query, set_list):
    for mol_set in set_list:
        mol_set.mol_list = mol_set.mol_list.filter(**build_molecules_filter(attribute, filter_query))



"""
this function is a helper function, building the actual filterfunction for filter_sets
"""


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
