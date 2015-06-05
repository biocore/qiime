"""Methods for computing microbial diversity indices

The basis for this library comes from Gevers et al. 2014
(http://www.ncbi.nlm.nih.gov/pubmed/24629344) in which a microbial dysbiosis
index was created based on observed increases and decreases in organisms with
respect to Crohn's disease.
"""
from __future__ import division

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2014, The QIIME project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


from numpy import log, nan, allclose, inf


def compute_index(table, increased, decreased, key):
    """Compute a per-sample index

    Parameters
    ----------
    table : biom.Table
        A biom table that has information associated with `key`, such as
        taxonomy.
    increased : set
        A set of items that have been observed to have increased
    decreased : set
        A set of items that have been observed to have decreased
    key : str
        The metadata key to use for the computation of the index.

    Returns
    -------
    generator
        (sample_id, index_score)

    Raises
    ------
    KeyError
        If the key isn't present
    ValueError
        If none of the increased or decreased items exist in the table

    Notes
    -----
    Yields `nan` if the decreased count is 0.

    """
    md_test = table.metadata(axis='observation')[0]
    if key not in md_test:
        raise KeyError("%s is not present, the following are present: %s" %
                       (key, ','.join(md_test)))

    inc_f = lambda v, i, md: set(md[key]) & increased
    dec_f = lambda v, i, md: set(md[key]) & decreased
    inc_t = table.filter(inc_f, axis='observation', inplace=False)
    dec_t = table.filter(dec_f, axis='observation', inplace=False)

    if inc_t.is_empty():
        raise ValueError("None of the increased items were found")

    if dec_t.is_empty():
        raise ValueError("None of the decreased items were found")

    ids_in_common = set(inc_t.ids()) & set(dec_t.ids())

    for id_ in ids_in_common:
        inc_count = inc_t.data(id_, dense=False).sum()
        dec_count = dec_t.data(id_, dense=False).sum()

        if allclose(dec_count, 0):
            yield (id_, nan)
        elif allclose(inc_count, 0):
            yield (id_, -inf)
        else:
            yield (id_, log(inc_count / dec_count))
