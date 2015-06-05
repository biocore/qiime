#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Michael Dwan", "Logan Knecht",
               "Damien Coy", "Levi McCracken"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import join
from types import ListType

import pandas as pd
from skbio.stats.distance import DistanceMatrix
from skbio.stats.distance import anosim, permanova, bioenv

from qiime.parse import parse_mapping_file_to_dict
from qiime.util import get_qiime_temp_dir, MetadataMap, RExecutor

methods = ['adonis', 'anosim', 'bioenv', 'morans_i', 'mrpp', 'permanova',
           'permdisp', 'dbrda']


def compare_categories(dm_fp, map_fp, method, categories, num_perms, out_dir):
    """Runs the specified statistical method using the category of interest.

    This method does not return anything; all output is written to results
    files in out_dir.

    Arguments:
        dm_fp - filepath to the input distance matrix
        map_fp - filepath to the input metadata mapping file
        categories - list of categories in the metadata mapping file to
            consider in the statistical test. Multiple categories will only be
            considered if method is 'bioenv', otherwise only the first category
            will be considered
        num_perms - the number of permutations to use when calculating the
            p-value. If method is 'bioenv' or 'morans_i', this parameter will
            be ignored as they are not permutation-based methods
        out_dir - path to the output directory where results files will be
            written. It is assumed that this directory already exists and we
            have write permissions to it
    """
    # Make sure we were passed a list of categories, not a single string.
    if not isinstance(categories, ListType):
        raise TypeError("The supplied categories must be a list of "
                        "strings.")

    # Special case: we do not allow SampleID as it is not a category, neither
    # in data structure representation nor in terms of a statistical test (no
    # groups are formed since all entries are unique IDs).
    if 'SampleID' in categories:
        raise ValueError("Cannot use SampleID as a category because it is a "
                         "unique identifier for each sample, and thus does "
                         "not create groups of samples (nor can it be used as "
                         "a numeric category in Moran's I or BIO-ENV "
                         "analyses). Please choose a different metadata "
                         "column to perform statistical tests on.")

    dm = DistanceMatrix.read(dm_fp)

    if method in ('anosim', 'permanova', 'bioenv'):
        with open(map_fp, 'U') as map_f:
            md_dict = parse_mapping_file_to_dict(map_f)[0]
        df = pd.DataFrame.from_dict(md_dict, orient='index')

        out_fp = join(out_dir, '%s_results.txt' % method)

        if method in ('anosim', 'permanova'):
            if method == 'anosim':
                method_fn = anosim
            elif method == 'permanova':
                method_fn = permanova

            results = method_fn(dm, df, column=categories[0],
                                permutations=num_perms)
        elif method == 'bioenv':
            results = bioenv(dm, df, columns=categories)

        results.to_csv(out_fp, sep='\t')
    else:
        # Remove any samples from the mapping file that aren't in the distance
        # matrix (important for validation checks). Use strict=True so that an
        # error is raised if the distance matrix contains any samples that
        # aren't in the mapping file.
        with open(map_fp, 'U') as map_f:
            md_map = MetadataMap.parseMetadataMap(map_f)
        md_map.filterSamples(dm.ids, strict=True)

        # These methods are run in R. Input validation must be done here before
        # running the R commands.
        if method in ['adonis', 'morans_i', 'mrpp', 'permdisp', 'dbrda']:
            # Check to make sure all categories passed in are in mapping file
            # and are not all the same value.
            for category in categories:
                if not category in md_map.CategoryNames:
                    raise ValueError("Category '%s' not found in mapping file "
                                     "columns." % category)

                if md_map.hasSingleCategoryValue(category):
                    raise ValueError("All values in category '%s' are the "
                                     "same. The statistical method '%s' "
                                     "cannot operate on a category that "
                                     "creates only a single group of samples "
                                     "(e.g. there are no 'between' distances "
                                     "because there is only a single group)."
                                     % (category, method))

            # Build the command arguments string.
            command_args = ['-d %s -m %s -c %s -o %s'
                            % (dm_fp, map_fp, categories[0], out_dir)]

            if method == 'morans_i':
                # Moran's I requires only numeric categories.
                for category in categories:
                    if not md_map.isNumericCategory(category):
                        raise TypeError("The category '%s' is not numeric. "
                                        "Not all values could be converted to "
                                        "numbers." % category)
            else:
                # The rest require groups of samples, so the category values
                # cannot all be unique.
                for category in categories:
                    if (md_map.hasUniqueCategoryValues(category) and not
                        (method == 'adonis' and
                         md_map.isNumericCategory(category))):
                        raise ValueError("All values in category '%s' are "
                                         "unique. This statistical method "
                                         "cannot operate on a category with "
                                         "unique values (e.g. there are no "
                                         "'within' distances because each "
                                         "group of samples contains only a "
                                         "single sample)." % category)

                # Only Moran's I doesn't accept a number of permutations.
                if num_perms < 0:
                    raise ValueError("The number of permutations must be "
                                     "greater than or equal to zero.")

                command_args[0] += ' -n %d' % num_perms

            rex = RExecutor(TmpDir=get_qiime_temp_dir())
            rex(command_args, '%s.r' % method)
        else:
            raise ValueError("Unrecognized method '%s'. Valid methods: %r"
                             % (method, methods))
