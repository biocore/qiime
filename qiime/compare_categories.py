#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Michael Dwan", "Logan Knecht",
               "Damien Coy", "Levi McCracken"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import join
from types import ListType

import pandas as pd
from skbio.core.distance import DistanceMatrix
from skbio.math.stats.distance import ANOSIM, PERMANOVA

from qiime.format import format_best_results
from qiime.parse import parse_mapping_file_to_dict
from qiime.stats import Best
from qiime.util import get_qiime_temp_dir, MetadataMap, RExecutor

methods = ['adonis', 'anosim', 'best', 'morans_i', 'mrpp', 'permanova',
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
            considered if method is 'best', otherwise only the first category
            will be considered
        num_perms - the number of permutations to use when calculating the
            p-value. If method is 'best' or 'morans_i', this parameter will be
            ignored as they are not permutation-based methods
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
                         "a numeric category in Moran's I or BEST analyses). "
                         "Please use a different metadata column to perform "
                         "statistical tests on.")

    with open(dm_fp, 'U') as dm_f:
        dm = DistanceMatrix.from_file(dm_f)

    # These methods are in skbio. There are still methods in qiime.stats that
    # need to be ported to skbio, at which point a lot of this logic can be
    # simplified.
    if method in ('anosim', 'permanova'):
        if method == 'anosim':
            method_cls = ANOSIM
        elif method == 'permanova':
            method_cls = PERMANOVA
        else:
            # Should never get here...
            pass

        with open(map_fp, 'U') as map_f:
            md_dict = parse_mapping_file_to_dict(map_f)[0]
        df = pd.DataFrame.from_dict(md_dict, orient='index')

        method_inst = method_cls(dm, df, column=categories[0])
        results = method_inst(num_perms)

        with open(join(out_dir, '%s_results.txt' % method), 'w') as out_f:
            out_f.write(results.summary())
    else:
        # Remove any samples from the mapping file that aren't in the distance
        # matrix (important for validation checks). Use strict=True so that an
        # error is raised if the distance matrix contains any samples that
        # aren't in the mapping file.
        with open(map_fp, 'U') as map_f:
            md_map = MetadataMap.parseMetadataMap(map_f)
        md_map.filterSamples(dm.ids, strict=True)

        # Run the specified statistical method.
        if method in ['adonis', 'morans_i', 'mrpp', 'permdisp', 'dbrda']:
            # These methods are run in R. Input validation must be done here
            # before running the R commands. The pure-Python implementations
            # perform all validation in the classes in the stats module.

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
                    if md_map.hasUniqueCategoryValues(category):
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
            rex(command_args, '%s.r' % method, output_dir=out_dir)
        elif method == 'best':
            best = Best(dm, md_map, categories)
            best_results = best()

            with open(join(out_dir, '%s_results.txt' % method), 'w') as out_f:
                out_f.write(format_best_results(best_results))
        else:
            raise ValueError("Unrecognized method '%s'. Valid methods: %r"
                             % (method, methods))
