"""Methods for computing microbial diversity indices

The basis for this library comes from Gevers et al. 2014
(http://www.ncbi.nlm.nih.gov/pubmed/24629344) in which a microbial dysbiosis
index was created based on observed increases and decreases in organisms with
respect to Crohn's disease.
"""

from __future__ import division

from numpy import log, nan


def compute_index(table, increased, decreased):
    """Compute a per-sample index

    Parameters
    ----------
    table : biom.Table
        A biom table that has taxonomy
    increased : set
        A set of taxon names that have been observed to have increased
    decreased : set
        A set of taxon names that have been observed to have decreased

    Raises
    ------
    AttributeError
        If taxonomy isn't present
    ValueError
        If none of the increased or decreased taxa exist in the table

    Notes
    -----
    Yields `nan` if both the decreased count is 0.

    Returns
    -------
    generator
        (sample_id, index_score)
    """
    if 'taxonomy' not in table.metadata(axis='observation')[0]:
        raise AttributeError("Missing taxonomy data")

    inc_f = lambda v, i, md: set(md['taxonomy']) & increased
    dec_f = lambda v, i, md: set(md['taxonomy']) & decreased
    inc_t = table.filter(inc_f, axis='observation', inplace=False)
    dec_t = table.filter(dec_f, axis='observation', inplace=False)

    if inc_t.is_empty():
        raise ValueError("None of the increased taxa were found")

    if dec_t.is_empty():
        raise ValueError("None of the decreased taxa were found")

    ids_in_common = set(inc_t.ids()) & (set(dec_t.ids()))

    for id_ in ids_in_common:
        inc_count = inc_t.data(id_, dense=False).sum()
        dec_count = dec_t.data(id_, dense=False).sum()

        if dec_count == 0:
            yield (id_, nan)
        else:
            yield (id_, log(inc_count / dec_count))
