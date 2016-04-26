#!/usr/bin/env python
# file make_rarefaction_plots.py
from __future__ import division
__author__ = "Meg Pirrung"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Meg Pirrung", "Jesse Stombaugh", "Antonio Gonzalez Pena",
               "Will Van Treuren", "Yoshiki Vazquez Baeza", "Jai Ram Rideout",
               "Evan Bolyen"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"

from matplotlib import use
use('Agg', warn=False)
from sys import exit
from qiime.parse import parse_rarefaction_data
from matplotlib.pyplot import savefig, clf, gca, gcf, errorbar
import matplotlib.pyplot as plt
import os.path
from os.path import splitext, split
from qiime.colors import iter_color_groups
from qiime.sort import natsort
from qiime.util import create_dir, stderr
from numpy import isnan, nan, array, transpose, mean, std, arange
from StringIO import StringIO
import urllib
import base64


def save_ave_rarefaction_plots(xaxis, yvals, err, xmax, ymax, ops,
                               mapping_category, imagetype, res, data_colors, colors, fpath,
                               background_color, label_color, metric_name, output_type="file_creation"):
    '''This function creates the images, using matplotlib.'''
    # Create the plot image
    plt.clf()
    plt.title(metric_name + ": " + mapping_category)
    fig = plt.gcf()

    # Add the lines to the plot
    for o in ops:

        l = o
        plt.errorbar(xaxis[:len(yvals[o])], yvals[o],
                     yerr=err[o][:len(yvals[o])], label=l, color=
                     data_colors[colors[o]].toHex(), elinewidth=1, lw=2, capsize=4)

    # get the plot axis
    ax = plt.gca()
    ax.set_axis_bgcolor(background_color)
    # ax.set_yscale('log',basey=2,basex=2)
    # set tick colors and width
    for line in ax.yaxis.get_ticklines():
        # line is a matplotlib.lines.Line2D instance
        line.set_color(label_color)
        line.set_markeredgewidth(1)

    for line in ax.xaxis.get_ticklines():
        # line is a matplotlib.lines.Line2D instance
        line.set_color(label_color)
        line.set_markeredgewidth(1)

    # set x/y limits and labels for plot
    ax.set_axisbelow(True)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.set_xlabel('Sequences Per Sample')
    ax.set_ylabel("Rarefaction Measure: " + metric_name)

    imgpath = fpath + mapping_category + '.' + imagetype
    if output_type == "file_creation":
        # Save the image
        plt.savefig(imgpath, format=imagetype, dpi=res)
        # Get the image name for the saved image relative to the main directory
        image_loc = imgpath
        plt.close()

        return
    elif (output_type == "memory"):
        imgdata = StringIO()
        plt.savefig(imgdata, format='png', dpi=res, transparent=True)
        imgdata.seek(0)
        plt.close()

        return {imgpath: imgdata}
    else:
        return None


def save_single_ave_rarefaction_plots(xaxis, yvals, err, xmax, ymax, ops,
                                      mapping_category, imagetype, res, data_colors, colors, fpath,
                                      background_color, label_color, rarefaction_legend_mat, metric_name,
                                      mapping_lookup, output_type="file_creation"):
    '''This function creates the images, using matplotlib.'''

    avg_plots = {}
    # Add the lines to the plot
    for o in ops:
        rarefaction_legend_mat[metric_name]['groups'][mapping_category][o]['ave_link'] = \
            os.path.join('html_plots',
                         metric_name + mapping_lookup[mapping_category + '-' + o] + '_ave.' + imagetype)

        # Create the plot image
        plt.clf()
        plt.title(metric_name + ": " + mapping_category, weight='regular')

        fig = plt.gcf()

        l = o
        plt.errorbar(xaxis[:len(yvals[o])], yvals[o],
                     yerr=err[o][:len(yvals[o])], label=l, color=
                     data_colors[colors[o]].toHex(), elinewidth=1, lw=2, capsize=4
                     )
        plt.alpha = (0)

        # get the plot axis
        ax = plt.gca()

        # set tick colors and width
        for line in ax.yaxis.get_ticklines():
            # line is a matplotlib.lines.Line2D instance
            line.set_color('black')
            line.set_markeredgewidth(1)

        for line in ax.xaxis.get_ticklines():
            # line is a matplotlib.lines.Line2D instance
            line.set_color('black')
            line.set_markeredgewidth(1)

        # set x/y limits and labels for plot
        ax.set_axisbelow(True)
        ax.set_xlim((0, xmax))
        ax.set_ylim((0, ymax))
        ax.set_xlabel('Sequences Per Sample')
        ax.set_ylabel("Rarefaction Measure: " + metric_name)

        x = ax.xaxis.get_label()
        x.set_weight('regular')
        # x.set_name('Arial')
        y = ax.yaxis.get_label()
        y.set_weight('regular')
        # y.set_name('Arial')

        imgpath = fpath + \
            mapping_lookup[mapping_category + '-' + o] + \
            '_ave.' + imagetype

        if output_type == "file_creation":
            # Save the image
            plt.savefig(
                imgpath,
                format=imagetype,
                dpi=res,
                transparent=True)
            # Get the image name for the saved image relative to the main
            # directory
            image_loc = imgpath
            plt.close()
        elif (output_type == "memory"):
            imgdata = StringIO()
            plt.savefig(imgdata, format='png', dpi=res, transparent=True)
            imgdata.seek(0)
            avg_plots[imgpath] = imgdata
            plt.close()

    if output_type == "file_creation":
        return rarefaction_legend_mat
    elif output_type == "memory":
        return rarefaction_legend_mat, avg_plots


def save_single_rarefaction_plots(sample_dict, imagetype, metric_name,
                                  data_colors, colors, fpath,
                                  background_color, label_color, res, ymax, xmax,
                                  rarefaction_legend_mat, groups,
                                  mapping_category, group_id, mapping_lookup, output_type="file_creation"):
    '''This function creates the images, using matplotlib.'''
    # Create the plot image
    plt.clf()
    # plt.title(str(metric_name))
    fig = plt.gcf()
    ax = fig.add_subplot(111)

    for o in groups:
        for i in sample_dict[o]:

            xaxis = []
            # this creates duplicates of the xval, since there are several
            # iterations
            for t in range(len(sample_dict[o][i])):
                xaxis.append(i)

            # If all the yvals are nan at a particular xval, skip adding
            # it to the plot
            if not isnan(sample_dict[o][i])[0]:
                scplot = ax.scatter(xaxis, sample_dict[o][i],
                                    c=data_colors[colors[o]].toHex(),
                                    marker='s', edgecolors='none')

    # get the plot axis
    ax = plt.gca()
    ax.set_axis_bgcolor(background_color)

    # set tick colors and width
    for line in ax.yaxis.get_ticklines():
        # line is a matplotlib.lines.Line2D instance
        line.set_color(label_color)
        line.set_markeredgewidth(1)

    for line in ax.xaxis.get_ticklines():
        # line is a matplotlib.lines.Line2D instance
        line.set_color(label_color)
        line.set_markeredgewidth(1)

    # set x/y limits and labels for plot
    ax.set_axisbelow(False)
    ax.set_xlim((0, xmax))
    ax.set_ylim((0, ymax))
    ax.set_xlabel('Sequences Per Sample')
    ax.set_ylabel("Rarefaction Measure: " + str(metric_name))

    x = ax.xaxis.get_label()
    x.set_weight('regular')
    # x.set_name('Arial')
    y = ax.yaxis.get_label()
    y.set_weight('regular')
    # y.set_name('Arial')

    # Create file for image
    imgpath = os.path.join(
        fpath,
        metric_name +
        mapping_lookup[
            mapping_category +
            '-' +
            group_id] +
        '_raw.' +
        imagetype)


    # Since both the average and raw are saved the same way we will save the
    # raw link as well
    rarefaction_legend_mat[metric_name]['groups'][mapping_category][group_id]['raw_link'] = \
        os.path.join('html_plots',
                     metric_name + mapping_lookup[mapping_category + '-' + group_id] + '_raw.' + imagetype)

    if output_type == "file_creation":
        # Save the image
        plt.savefig(imgpath, format=imagetype, dpi=res, transparent=True)
        # Get the image name for the saved image relative to the main directory
        image_loc = imgpath
        plt.close()

        return rarefaction_legend_mat
    elif (output_type == "memory"):
        imgdata = StringIO()
        plt.savefig(imgdata, format='png', dpi=res, transparent=True)
        imgdata.seek(0)
        plt.close()

        return [rarefaction_legend_mat, {imgpath: imgdata}]


def get_rarefaction_data(rarefaction_data, col_headers):
    '''This function takes a rarefaction file and converts it into an array'''

    rare_mat_raw = array(rarefaction_data)
    rare_mat_min = [rare_mat_raw[x][2:] for x in range(0, len(rare_mat_raw))]
    seqs_per_samp = [rare_mat_raw[x][0] for x in range(0, len(rare_mat_raw))]
    sampleIDs = col_headers[3:]

    # Need to transpose the array to be used in averaging
    rare_mat_trans = transpose(array(rare_mat_min)).tolist()

    return rare_mat_trans, seqs_per_samp, sampleIDs


def ave_seqs_per_sample(matrix, seqs_per_samp, sampleIDs):
    """Calculate the average for each sampleID across each number of \
    seqs/sample"""

    ave_ser = {}
    temp_dict = {}
    # Iterate through the samples id's and create a dictionary
    for i, sid in enumerate(sampleIDs):
        temp_dict[sid] = {}
        for j, seq in enumerate(seqs_per_samp):
            try:
                temp_dict[sid][seq].append(matrix[i][j])
            except(KeyError):
                temp_dict[sid][seq] = []
                temp_dict[sid][seq].append(matrix[i][j])

    # create a dictionary for average data
    for sid in sampleIDs:
        ave_ser[sid] = []
        keys = sorted(temp_dict[sid].keys())

        for k in keys:
            ave_ser[sid].append(mean(array(temp_dict[sid][k]), 0))

    return ave_ser


def make_error_series(rare_mat, groups, std_type):
    """Create mean and error bar series for the supplied mapping category"""

    err_ser = dict()
    collapsed_ser = dict()
    seen = set()
    pre_err = {}

    ops = [k for k in groups]

    notfound = []
    # Iterate through the groups
    for o in ops:
        pre_err[o] = []

        # For each sample in group, create a row in a list
        for samID in groups[o]:
            pre_err[o].append(rare_mat[samID])

        min_len = min([len(i) - i.count('nan') for i in pre_err[o]])
        pre_err[o] = [x[:min_len] for x in pre_err[o]]

    # iterate through the groups and calculate std deviations and error
    for o in ops:
        opsarray = array(pre_err[o])
        mn = mean(opsarray, 0)
        collapsed_ser[o] = mn.tolist()

        if std_type == 'stderr':
            # this calculates the standard error
            # (using sample standard deviation)
            stderr_result = stderr(opsarray, 0)
            err_ser[o] = stderr_result.tolist()
        else:
            # this calculates the population standard deviation
            stddev = std(opsarray, 0)
            err_ser[o] = stddev.tolist()

    return collapsed_ser, err_ser, ops


def save_rarefaction_data(rare_mat, xaxis, xmax,
                          mapping_category, colors, rare_type, data_colors, groups, std_type):
    '''This function formats the average data and writes it to the output
       directory'''

    # get the error data
    yaxis, err, ops = make_error_series(rare_mat, groups, std_type)

    lines = []
    lines.append("# " + rare_type + '\n')
    lines.append("# " + mapping_category + '\n')
    line = ''
    line += 'xaxis: '
    for v in xaxis:
        line += str(v) + '\t'
    line += '\n'
    lines.append(line)
    lines.append('xmax: ' + str(xmax) + '\n')

    for o in colors.keys():
        lines.append(">> " + o + '\n')

        # write the color lines
        if colors is not None:
            try:
                lines.append("color " + data_colors[colors[o]].toHex() + '\n')
            except(KeyError):
                print 'Color reference is missing!'

        # write the rarefection series lines
        lines.append('series ')
        line = ''
        try:
            for v in yaxis[o]:
                line += str(v) + '\t'
        except(TypeError):
            line += str(yaxis[o])
        line += '\n'
        lines.append(line)

        # write the rarefaction error lines
        lines.append('error ')
        line = ''
        try:
            for e in err[o]:
                if e == 0:
                    line += str(nan) + '\t'
                else:
                    line += str(e) + '\t'

        except(TypeError):
            line += str(err[o])
        line += '\n'
        lines.append(line)

    return lines


def make_averages(color_prefs, data, background_color, label_color, rares,
                  output_dir, resolution, imagetype, ymax, suppress_webpage,
                  std_type, output_type="file_creation",
                  generate_per_sample_plots=True,
                  generate_average_tables=True):
    '''This is the main function, which takes the rarefaction files, calls the
        functions to make plots and formatting the output html.'''
    rarelines = []
    rarefaction_legend_mat = {}

    if ymax:
        user_ymax = True
    else:
        user_ymax = False

    if not suppress_webpage and output_type == "file_creation":
        # in this option the path must include the output directory
        ave_output_dir = os.path.join(output_dir, 'average_plots')
        if generate_per_sample_plots:
            all_output_dir = os.path.join(output_dir, 'html_plots')
            # Create the directories, where plots and data will be written
            create_dir(all_output_dir)
        else:
            all_output_dir = ""

    elif output_type == 'memory':
        # this is rather an artificial path to work with the javascript code
        ave_output_dir = 'plot/average_plots'
        if generate_per_sample_plots:
            all_output_dir = 'plot/html_plots'
        else:
            all_output_dir = ""

    ave_data_file_path = os.path.join(output_dir, 'average_tables')
    if output_type == "file_creation":
        create_dir(ave_output_dir)
        if generate_average_tables:
            create_dir(ave_data_file_path, False)

    metric_num = 0
    rarefaction_legend_mat = {}
    rarefaction_data_mat = {}
    rare_num = 0

    # this is  a fix for the issue of writing field values as the filenames
    mapping_lookup = {}
    for i, column in enumerate(data['map'][0]):
        for j, row in enumerate(data['map'][1:]):
            mapping_lookup['%s-%s' % (column, row[i])] = 'col_%s_row_%s' % \
                (str(i), str(j))

    all_plots = []
    # Iterate through the rarefaction files
    for r in natsort(rares):
        raredata = rares[r]
        metric_name = r.split('.')[0]

        # convert the rarefaction data into variables
        col_headers, comments, rarefaction_fn, rarefaction_data = rares[r]

        # Here we only need to perform these steps once, since the data is
        # the same for all rarefaction files
        if rare_num == 0:

            # Remove samples from the mapping file, which contain no data after
            # rarefaction
            updated_mapping = []
            for j in data['map']:

                # Add the mapping header
                if j[0] == 'SampleID':
                    updated_mapping.append(j)

                # Determine if the sample exists in the rarefaction file
                for i in col_headers[3:]:
                    if j[0] == i:
                        updated_mapping.append(j)

            # Get the groups and colors for the updated mapping file
            groups_and_colors = iter_color_groups(updated_mapping, color_prefs)
            groups_and_colors = list(groups_and_colors)

        # parse the rarefaction data

        rare_mat_trans, seqs_per_samp, sampleIDs = \
            get_rarefaction_data(rarefaction_data, col_headers)

        rarefaction_legend_mat[metric_name] = {}

        # Create dictionary variables and get the colors for each Sample
        sample_colors = None
        rarefaction_legend_mat[metric_name]['groups'] = {}
        for i in range(len(groups_and_colors)):
            labelname = groups_and_colors[i][0]
            # Create a legend dictionary for html output
            rarefaction_legend_mat[metric_name]['groups'][labelname] = {}
            # If this is the first time iterating through the rarefaction data
            # create a data dictionary for html output
            if rare_num == 0:
                rarefaction_data_mat[labelname] = {}

            # If the labelname is SampleID, use the colors assigned
            if labelname == 'SampleID':
                sample_colors = groups_and_colors[i][2]
                sample_data_colors = groups_and_colors[i][3]

        rare_num = 1

        # If sample colors were not assigned, create a list of sample colors
        if not sample_colors:
            samples_and_colors = iter_color_groups(updated_mapping,
                                                   {'SampleID': {'column': 'SampleID', 'colors':
                                                                 (('red', (0, 100, 100)), ('blue', (240, 100, 100)))}})
            samples_and_colors = list(samples_and_colors)
            sample_colors = samples_and_colors[0][2]
            sample_data_colors = samples_and_colors[0][3]

        sample_dict = {}
        # Create a dictionary containing the samples
        for i, sid in enumerate(sampleIDs):
            if sid in (i[0] for i in updated_mapping):
                sample_dict[sid] = {}
                for j, seq in enumerate(seqs_per_samp):
                    try:
                        sample_dict[sid][seq].append(rare_mat_trans[i][j])
                    except(KeyError):
                        sample_dict[sid][seq] = []
                        sample_dict[sid][seq].append(rare_mat_trans[i][j])

        # convert xvals to float
        xaxisvals = sorted([float(x) for x in set(seqs_per_samp)])

        # get the rarefaction averages
        rare_mat_ave = ave_seqs_per_sample(rare_mat_trans, seqs_per_samp,
                                           sampleIDs)

        # calculate the max xval
        xmax = max(xaxisvals) + (xaxisvals[len(xaxisvals) - 1] -
                                 xaxisvals[len(xaxisvals) - 2])

        if not user_ymax:
            ymax = 0
            for i in range(len(groups_and_colors)):
                labelname = groups_and_colors[i][0]
                groups = groups_and_colors[i][1]
                colors = groups_and_colors[i][2]
                data_colors = groups_and_colors[i][3]
                if generate_average_tables:
                    ave_file_path = os.path.join(
                        ave_data_file_path,
                        metric_name)
                    # save the rarefaction averages

                rare_lines = save_rarefaction_data(
                    rare_mat_ave, xaxisvals, xmax,
                    labelname, colors, r, data_colors, groups,
                    std_type)

                # write out the rarefaction average data
                if output_type == "file_creation" and generate_average_tables:
                    open(
                        ave_file_path +
                        labelname +
                        '.txt',
                        'w').writelines(
                        rare_lines)

                # take the formatted rarefaction averages and format the
                # results
                rares_data = parse_rarefaction_data(
                    ''.join(rare_lines[:]).split('\n'))

                # determine the ymax based on the average data
                # multiple the ymax, since the dots can end up on the border
                new_ymax = (max([max(v) for v in rares_data['series'].values()]) +
                            max([max(e) for e in rares_data['error'].values()])) * 1.15
                if isnan(new_ymax):
                    new_ymax = (max([max(v) for v in
                                    rares_data['series'].values()])) * 1.15

                if new_ymax > ymax:
                    ymax = new_ymax

        iterator_num = 0

        # iterate through the groups
        for i in range(len(groups_and_colors)):
            labelname = groups_and_colors[i][0]
            groups = groups_and_colors[i][1]
            colors = groups_and_colors[i][2]
            data_colors = groups_and_colors[i][3]
            data_color_order = groups_and_colors[i][4]

            # save the rarefaction averages
            rare_lines = save_rarefaction_data(rare_mat_ave, xaxisvals, xmax,
                                               labelname, colors, r, data_colors, groups,
                                               std_type)

            # take the formatted rarefaction averages and format the results
            rares_data = parse_rarefaction_data(
                ''.join(rare_lines[:]).split('\n'))

            if not suppress_webpage:

                if iterator_num == 0:
                    rarefaction_legend_mat[metric_name]['samples'] = {}
                    for o in sample_dict:
                        rarefaction_legend_mat[metric_name]['samples'][o] = {}
                        # Add values to the legend dictionary
                        rarefaction_legend_mat[metric_name]['samples'][o][
                            'color'] = sample_data_colors[sample_colors[o]].toHex()

                    iterator_num = 1

                # Iterate through the groups and create the legend dictionary
                for g in groups:
                    # create a dictionary of samples and their colors
                    rarefaction_legend_mat[
                        metric_name][
                        'groups'][
                        labelname][
                        g] = {
                    }
                    rarefaction_legend_mat[metric_name]['groups'][
                        labelname][g]['groupsamples'] = groups[g]
                    rarefaction_legend_mat[metric_name]['groups'][labelname][g]['groupcolor'] =\
                        data_colors[colors[g]].toHex()

                # Create the individual category average plots
                if output_type == "file_creation":
                    rarefaction_data_mat, rarefaction_legend_mat = make_plots(
                        background_color, label_color,
                        rares_data, ymax, xmax, all_output_dir,
                        resolution, imagetype, groups, colors,
                        data_colors, metric_name, labelname,
                        rarefaction_data_mat, rarefaction_legend_mat,
                        sample_dict, sample_data_colors,
                        sample_colors, mapping_lookup, output_type,
                        generate_per_sample_plots)
                elif output_type == "memory":
                    rarefaction_data_mat, rarefaction_legend_mat, all_plots_single, \
                        all_plots_ave = make_plots(
                            background_color, label_color,
                            rares_data, ymax, xmax, all_output_dir,
                            resolution, imagetype, groups, colors,
                            data_colors, metric_name, labelname,
                            rarefaction_data_mat, rarefaction_legend_mat,
                            sample_dict, sample_data_colors,
                            sample_colors, mapping_lookup, output_type,
                            generate_per_sample_plots)

                # generate the filepath for the image file
                file_path = os.path.join(ave_output_dir,
                                         splitext(split(rares_data['headers'][0])[1])[0])

                # Create the average plots
                categories = [k for k in groups]
                all_plots_rare = save_ave_rarefaction_plots(
                    rares_data['xaxis'], rares_data['series'],
                    rares_data[
                        'error'], xmax, ymax, categories,
                    labelname, imagetype, resolution, data_colors,
                    colors, file_path, background_color, label_color,
                    metric_name, output_type)

                if output_type == "memory":
                    all_plots.append(all_plots_rare)
                    all_plots.extend(all_plots_single)
                    all_plots.append(all_plots_ave)
            else:
                # generate the filepath for the image file
                file_path = os.path.join(ave_output_dir,
                                         splitext(split(rares_data['headers'][0])[1])[0])

                categories = [k for k in groups]
                all_plots_rare = save_ave_rarefaction_plots(
                    rares_data['xaxis'], rares_data['series'],
                    rares_data[
                        'error'], xmax, ymax, categories,
                    labelname, imagetype, resolution, data_colors,
                    colors, file_path, background_color, label_color,
                    metric_name, output_type)

    if not suppress_webpage:
        # format the html output
        html_output = make_html(rarefaction_legend_mat,
                                rarefaction_data_mat, xaxisvals, imagetype, mapping_lookup,
                                output_type, all_plots, generate_per_sample_plots)
    else:
        html_output = None

    return html_output


def make_html(rarefaction_legend_mat, rarefaction_data_mat, xaxisvals,
              imagetype, mapping_lookup, output_type="file_creation", all_plots=None,
              generate_per_sample_plots=True):
    rarefaction_legend_mat
    legend_td = [
        '<b>Legend</b><div STYLE="border: thin black solid; height: 300px; width: 200px; font-size: 12px; overflow: auto;"><table>']
    summarized_table = []
    metric_select_html = []
    category_select_html = []
    data_table_html = []
    metrics = []
    category_colors = {}
    cat_iter = 0
    # iterate the legend dictionary
    for m in natsort(rarefaction_legend_mat):

        # Create the metric select box options
        metric_select_html.append('<option value="%s">%s</option>' % (m, m))
        metrics.append(m)

        # iterate through the categories in the legend dictionary
        for category in natsort(rarefaction_legend_mat[m]['groups']):
            category_colors[category] = {}
            # Create the select box options
            if cat_iter == 0:
                cat_links = []
                for i in rarefaction_legend_mat[m]['groups'][category]:
                    cat_links.append(mapping_lookup[category + '-' + i])
                category_select_html.append('<option value="%s">%s</option>' %
                                            (category + '$#!' + '$#!'.join(cat_links), category))

            plot_iterator = 0
            # iterate through the groups in the legend dictionary and create
            # the html formatted rows for each category and group
            for group in natsort(rarefaction_legend_mat[m]['groups'][category]):
                sample_list = []
                category_colors[category][group] =\
                    rarefaction_legend_mat[m]['groups'][
                        category][group]['groupcolor']

                for sample in natsort(rarefaction_legend_mat[m]['groups'][category][group]['groupsamples']):
                    sample_list.append('\'' + sample + '\'')

                plot_iterator = plot_iterator + 1

                if generate_per_sample_plots:
                    legend_td.append(
                        '<tr id="%s" name="%s" style="display: none;"><td class="data" onmouseover="document.body.style.cursor=\'pointer\'"  onmouseout="document.body.style.cursor=\'default\'" onclick="toggle(%s)" id="%s" name="%s">&#x25B6;</td><td><input name="%s" type="checkbox" checked="True" onclick="show_hide_category(this)"></td><td style="color:%s">&#x25A0;&nbsp;</td><td class="data"><b>%s</b></td></tr>' % (m + category,
                                                                                                                                                                                                                                                                                                                                                                                                                                  m +
                                                                                                                                                                                                                                                                                                                                                                                                                                  category,
                                                                                                                                                                                                                                                                                                                                                                                                                                  "'" + m +
                                                                                                                                                                                                                                                                                                                                                                                                                                  mapping_lookup[
                                                                                                                                                                                                                                                                                                                                                                                                                                      category + '-' + group] + "'",
                                                                                                                                                                                                                                                                                                                                                                                                                                  m +
                                                                                                                                                                                                                                                                                                                                                                                                                                  mapping_lookup[
                                                                                                                                                                                                                                                                                                                                                                                                                                      category +
                                                                                                                                                                                                                                                                                                                                                                                                                                      '-' +
                                                                                                                                                                                                                                                                                                                                                                                                                                      group],
                                                                                                                                                                                                                                                                                                                                                                                                                                  ','.join(
                                                                                                                                                                                                                                                                                                                                                                                                                                      sample_list),
                                                                                                                                                                                                                                                                                                                                                                                                                                  m + mapping_lookup[category + '-' + group] +
                                                                                                                                                                                                                                                                                                                                                                                                                                  '_raw.' +
                                                                                                                                                                                                                                                                                                                                                                                                                                  imagetype,
                                                                                                                                                                                                                                                                                                                                                                                                                                  rarefaction_legend_mat[m]['groups'][
                                                                                                                                                                                                                                                                                                                                                                                                                                      category][
                                                                                                                                                                                                                                                                                                                                                                                                                                      group][
                                                                                                                                                                                                                                                                                                                                                                                                                                      'groupcolor'],
                                                                                                                                                                                                                                                                                                                                                                                                                                  group))
                else:
                    legend_td.append(
                        '<tr id="%s" name="%s" style="display: none;"><td class="data" onmouseover="document.body.style.cursor=\'pointer\'"  onmouseout="document.body.style.cursor=\'default\'" onclick="toggle(%s)" id="%s" name="%s">&#x25B6;</td><td>&nbsp;</td><td style="color:%s">&#x25A0;&nbsp;</td><td class="data"><b>%s</b></td></tr>' % (m + category,
                                                                                                                                                                                                                                                                                                                                                     m +
                                                                                                                                                                                                                                                                                                                                                     category,
                                                                                                                                                                                                                                                                                                                                                     "'" + m +
                                                                                                                                                                                                                                                                                                                                                     mapping_lookup[
                                                                                                                                                                                                                                                                                                                                                         category + '-' + group] + "'",
                                                                                                                                                                                                                                                                                                                                                     m +
                                                                                                                                                                                                                                                                                                                                                     mapping_lookup[
                                                                                                                                                                                                                                                                                                                                                         category +
                                                                                                                                                                                                                                                                                                                                                         '-' +
                                                                                                                                                                                                                                                                                                                                                         group],
                                                                                                                                                                                                                                                                                                                                                     ','.join(
                                                                                                                                                                                                                                                                                                                                                         sample_list),
                                                                                                                                                                                                                                                                                                                                                     rarefaction_legend_mat[m]['groups'][
                                                                                                                                                                                                                                                                                                                                                         category][
                                                                                                                                                                                                                                                                                                                                                         group][
                                                                                                                                                                                                                                                                                                                                                         'groupcolor'],
                                                                                                                                                                                                                                                                                                                                                     group))

                for sample in natsort(rarefaction_legend_mat[m]['groups'][category][group]['groupsamples']):
                    sample = str(sample)
                    legend_td.append(
                        '<tr id="%s" name="%s" style="display: none;"><td class="data" align="right">&#x221F;</td><td></td><td style="color:%s">&#x25C6;</td><td class="data" align="left"><b>%s</b></td></tr>' %
                        (m + mapping_lookup[category + '-' + group] + '_raw', m + mapping_lookup[category + '-' + group], rarefaction_legend_mat[m]['samples'][sample]['color'], sample))

        cat_iter = 1

    # iterate through the data dictionary and format the rows for the html
    # data table
    for category in rarefaction_data_mat:
        data_table_html.append(
            '<tr name="%s" style="display: none;"><td class="headers">%s</td><td class="headers">Seqs/Sample</td>' %
            (category, category))
        for j in metrics:
            data_table_html.append(
                '<td class="headers">%s Ave.</td><td class="headers">%s Err.</td>' %
                (j, j))
        data_table_html.append('</tr>')
        #data_table_html.append('<tr name="%s" style="display: none;"></tr>' % (category))
        for g in natsort(rarefaction_data_mat[category]):
            for i in range(len(xaxisvals)):
                data_table_html.append(
                    '<tr name="%s" style="display: none;">' %
                    (category))
                data_table_html.append(
                    '<td class="data" bgcolor="%s">%s</td><td class="data">%s</td>' %
                    (category_colors[category][g], g, xaxisvals[i]))
                # bugfix, was rarefaction_data_mat[category][g]
                for m in metrics:
                    data_table_html.append(
                        '<td class="data">%s</td><td class="data">%s</td>' %
                        (rarefaction_data_mat[category][g][m]['ave'][i], rarefaction_data_mat[category][g][m]['err'][i]))
        data_table_html.append('</tr>')

    legend_td.append('</table></div></div>')
    # Create the table that contains the plots and table
    plot_html = '%s' % ('\n'.join(legend_td))

    if output_type == "file_creation":
        # insert the formatted rows into the html string at the bottom of this
        # file
        if generate_per_sample_plots:
            html_output = HTML % ('',
                                  "img.setAttribute('src',\"./html_plots/\"+SelObject.value+array[i]+'_ave'+imagetype)",
                                  "img.setAttribute('src',\"./html_plots/\"+metric+array[i]+'_ave'+imagetype)",
                                  "img.setAttribute('src',\"./html_plots/\"+arguments[0]+'_raw'+imagetype)",
                                  '.' + imagetype,
                                  '\n'.join(metric_select_html),
                                  '\n'.join(category_select_html),
                                  plot_html,
                                  '\n'.join(data_table_html))
        else:
            html_output = HTML % ('',
                                  "img.setAttribute('src',\"./average_plots/\"+SelObject.value+array[0]+imagetype)",
                                  "img.setAttribute('src',\"./average_plots/\"+metric+array[0]+imagetype)",
                                  "",
                                  '.' + imagetype,
                                  '\n'.join(metric_select_html),
                                  '\n'.join(category_select_html),
                                  plot_html,
                                  '\n'.join(data_table_html))
    elif output_type == "memory":
        plots_html = ['all_plots = {}']
        for elements in all_plots:
            for k, v in elements.items():
                # the path is compatible with the javascript, see make_averages
                plots_html.append('all_plots["%s"] = "%s"' % (k,
                                                              "data:image/png;base64," + urllib.quote(base64.b64encode(v.buf))))

        # insert the formatted rows into the html string at the bottom of this
        # file
        if generate_per_sample_plots:
            html_output = HTML % ('\n'.join(plots_html),
                                  "img.setAttribute('src',all_plots[\"plot/html_plots/\"+SelObject.value+array[i]+'_ave'+imagetype])",
                                  "img.setAttribute('src',all_plots[\"plot/html_plots/\"+metric+array[i]+'_ave'+imagetype])",
                                  "img.setAttribute('src',all_plots[\"plot/html_plots/\"+arguments[0]+'_raw'+imagetype])",
                                  '.' + imagetype,
                                  '\n'.join(metric_select_html),
                                  '\n'.join(category_select_html),
                                  plot_html,
                                  '\n'.join(data_table_html))
        else:
            html_output = HTML % ('',
                                  "img.setAttribute('src',\"./average_plots/\"+SelObject.value+array[0]+imagetype)",
                                  "img.setAttribute('src',\"./average_plots/\"+metric+array[0]+imagetype)",
                                  "",
                                  '.' + imagetype,
                                  '\n'.join(metric_select_html),
                                  '\n'.join(category_select_html),
                                  plot_html,
                                  '\n'.join(data_table_html))
    return html_output


def make_plots(background_color, label_color, rares, ymax, xmax,
               output_dir, resolution, imagetype, groups, colors, data_colors,
               metric_name, labelname, rarefaction_data_mat,
               rarefaction_legend_mat, sample_dict, sample_data_colors,
               sample_colors, mapping_lookup, output_type="file_creation",
               generate_per_sample_plots=True):
    '''This is the main function for generating the rarefaction plots and html
        file.'''
    # Get the alpha rare data
    raredata = rares

    # generate the filepath for the image file
    file_path = os.path.join(output_dir,
                             splitext(split(raredata['headers'][0])[1])[0])

    all_plots_single = []
    # Sort and iterate through the groups
    for i in natsort(groups):
        # for k in groups[i]:
        for j in range(len(raredata['xaxis'])):
            group_field = i

            seq_per_sample_field = int(raredata['xaxis'][j])
            color_field = data_colors[colors[group_field]].toHex()

            # If a field is missing, then it means that one of the
            # samples did not contain enough sequences.
            # For this case, we will assign the value as n.a.
            try:
                average_field = raredata['series'][i][j]
                error_field = raredata['error'][i][j]
                if isnan(average_field):
                    error_field = nan
            except:
                average_field = nan
                error_field = nan

            # Add context to the data dictionary, which will be used in the
            # html
            if i in rarefaction_data_mat[labelname]:
                if metric_name in rarefaction_data_mat[labelname][i]:
                    rarefaction_data_mat[labelname][i][metric_name]['ave'].append(
                        ''.join('%10.3f' %
                                ((raredata['series'][i][j]))))
                    rarefaction_data_mat[labelname][i][metric_name]['err'].append(
                        ''.join('%10.3f' %
                                ((raredata['error'][i][j]))))
                else:
                    rarefaction_data_mat[labelname][i][metric_name] = {}
                    rarefaction_data_mat[labelname][i][metric_name]['ave'] = []
                    rarefaction_data_mat[labelname][i][metric_name]['err'] = []
                    rarefaction_data_mat[labelname][i][metric_name]['ave'].append(
                        ''.join('%10.3f' %
                                ((raredata['series'][i][j]))))
                    rarefaction_data_mat[labelname][i][metric_name]['err'].append(
                        ''.join('%10.3f' %
                                ((raredata['error'][i][j]))))

            else:
                rarefaction_data_mat[labelname][i] = {}
                rarefaction_data_mat[labelname][i][metric_name] = {}
                rarefaction_data_mat[labelname][i][metric_name]['ave'] = []
                rarefaction_data_mat[labelname][i][metric_name]['err'] = []
                rarefaction_data_mat[labelname][i][metric_name]['ave'].append(
                    ''.join('%10.3f' %
                            ((raredata['series'][i][j]))))
                rarefaction_data_mat[labelname][i][metric_name]['err'].append(
                    ''.join('%10.3f' %
                            ((raredata['error'][i][j]))))

        # Create raw plots for each group in a category
        if generate_per_sample_plots:
            if output_type == "file_creation":
                rarefaction_legend_mat = save_single_rarefaction_plots(
                    sample_dict,
                    imagetype, metric_name,
                    sample_data_colors, sample_colors,
                    output_dir, background_color,
                    label_color, resolution, ymax, xmax,
                    rarefaction_legend_mat, groups[i],
                    labelname, i, mapping_lookup, output_type)
            elif output_type == "memory":
                rarefaction_legend_mat, rare_plot_for_all = save_single_rarefaction_plots(
                    sample_dict,
                    imagetype, metric_name,
                    sample_data_colors, sample_colors,
                    output_dir, background_color,
                    label_color, resolution, ymax, xmax,
                    rarefaction_legend_mat, groups[i],
                    labelname, i, mapping_lookup, output_type)
                all_plots_single.append(rare_plot_for_all)

    all_plots_ave = {}
    if generate_per_sample_plots:
        # Create the rarefaction average plot and get updated legend information
        categories = [k for k in groups]

        if output_type == "file_creation":
            rarefaction_legend_mat = save_single_ave_rarefaction_plots(
                raredata['xaxis'],
                raredata['series'], raredata[
                    'error'], xmax, ymax, categories,
                labelname, imagetype, resolution, data_colors,
                colors, file_path, background_color, label_color,
                rarefaction_legend_mat, metric_name, mapping_lookup, output_type)
        elif output_type == "memory":
            rarefaction_legend_mat, all_plots_ave = save_single_ave_rarefaction_plots(
                raredata['xaxis'],
                raredata['series'], raredata[
                    'error'], xmax, ymax, categories,
                labelname, imagetype, resolution, data_colors,
                colors, file_path, background_color, label_color,
                rarefaction_legend_mat, metric_name, mapping_lookup, output_type)

    if output_type == "file_creation":
        return rarefaction_data_mat, rarefaction_legend_mat
    elif output_type == "memory":
        return (
            rarefaction_data_mat, rarefaction_legend_mat, all_plots_single, all_plots_ave
        )


HTML = '''
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
  <meta http-equiv="content-type" content="text/html;">
  <title>Rarefaction Curves</title>
<style type="text/css">
td.data{font-size:10px;border-spacing:0px 10px;text-align:center;}
td.headers{font-size:12px;font-weight:bold;text-align:center;}
table{border-spacing:0px;}
.removed{display:none;}
.expands{cursor:pointer; cursor:hand;}
.child1 td:first-child{padding-left: 3px;}
</style>
<script language="javascript" type="text/javascript">
%s
function show_hide_category(checkobject){
    var imagetype=document.getElementById('imagetype').value;
    img=document.getElementById(checkobject.name.replace('_raw'+imagetype,'_ave'+imagetype))
    if (checkobject.checked==false){
        img.style.display='none';
    }else{
        img.style.display='';
    }
}

function reset_tree(){
    var category=document.getElementById('category').value;
    var metric=document.getElementById('metric').value;
    var old_all_categories=document.getElementById('all_categories');
    var imagetype=document.getElementById('imagetype').value;
    cat_list=old_all_categories.value.split('$#!')
    if (metric!='' && category != ''){
    for (var i=1, il=cat_list.length; i<il; i++){
        group=metric+category+cat_list[i]
        main_class=metric+category
        var exp_item=document.getElementById(group);
        if (exp_item!=null){
            if (exp_item.innerHTML=='\u25BC'){
                exp_item.innerHTML='\u25B6'
                var rows=document.getElementsByName(group);
                for (var j=0, jl=rows.length; j<jl; j++){
                    rows[j].style.display="none";
                }
            }
            var rows=document.getElementsByName(group+'_raw'+imagetype);
            for (var j=0, jl=rows.length; j<jl; j++){
                if (rows[j].checked==false){
                    rows[j].checked=true;
                }
            }
        }
    }
}
}

function changeMetric(SelObject){
    var category=document.getElementById('category');
    var old_metric=document.getElementById('metric');
    var imagetype=document.getElementById('imagetype').value;
    var legend=document.getElementById('legend');
    var array=document.getElementById('all_categories').value.split('$#!')
    var plots=document.getElementById('plots');
    plots.style.display='none'
    reset_tree();
    if (category.value != ''){
        legend.style.display="";
        cat=SelObject.value+category.value
        data_display=document.getElementsByName(cat)
        for (var i=0, il=data_display.length; i<il; i++){
            data_display[i].style.display="";
        }
        cat=old_metric.value+category.value
        data_hide=document.getElementsByName(cat)
        for (var i=0, il=data_hide.length; i<il; i++){
            data_hide[i].style.display="none";
        }
        data_display=document.getElementsByName(category.value)
        for (var i=0, il=data_display.length; i<il; i++){
            data_display[i].style.display="";
        }
        new_cat=SelObject.value+category.value
        plots.innerHTML=''
        for (var i=1, il=array.length; i<il; i++){
            img=document.createElement('img')
            img.setAttribute('width',"600px")
            img.setAttribute('id',array[i]+'_ave'+imagetype)
            img.setAttribute('style','position:absolute;z-index:0')
            %s
            plots.appendChild(img)
        }
        plots.style.display=''
    }

old_metric.value=SelObject.value;

// If both combo boxes have changed the value, display a disclaimer
if (document.getElementById('select_metric_combo').selectedIndex !== 0 &&  document.getElementById('select_category_combo').selectedIndex !== 0) {
    document.getElementById('nan_disclaimer').style.display='inline';
}
}

function changeCategory(SelObject){
    var old_category=document.getElementById('category');
    var metric=document.getElementById('metric').value;
    var imagetype=document.getElementById('imagetype').value;
    var legend=document.getElementById('legend');
    var plots=document.getElementById('plots');
    var array=SelObject.value.split('$#!')
    var old_all_categories=document.getElementById('all_categories')
    category=array[0]
    plots.style.display='none'
    reset_tree();

    if (metric != ''){
        legend.style.display="";

        data_display=document.getElementsByName(category)
        for (var i=0, il=data_display.length; i<il; i++){
            data_display[i].style.display="";
        }
        data_hide=document.getElementsByName(old_category.value)
        for (var i=0, il=data_hide.length; i<il; i++){
            data_hide[i].style.display="none";
        }
        cat=metric+category
        data_display=document.getElementsByName(cat)
        for (var i=0, il=data_display.length; i<il; i++){
            data_display[i].style.display="";
        }
        cat=metric+old_category.value
        data_hide=document.getElementsByName(cat)
        for (var i=0, il=data_hide.length; i<il; i++){
            data_hide[i].style.display="none";
        }
        cat=metric+category
        plots.innerHTML=''
        for (var i=1, il=array.length; i<il; i++){
            img=document.createElement('img')
            img.setAttribute('width',"600px")
            img.setAttribute('id',metric+array[i]+'_ave'+imagetype)
            img.setAttribute('style','position:absolute;z-index:0')
            %s
            plots.appendChild(img)
        }
        plots.style.display=''
    }
old_all_categories.value=SelObject.value;
old_category.value=category;
// If both combo boxes have changed the value, display a disclaimer
if (document.getElementById('select_metric_combo').selectedIndex !== 0 &&  document.getElementById('select_category_combo').selectedIndex !== 0) {
    document.getElementById('nan_disclaimer').style.display='inline';
}
}
function toggle(){
    var plots=document.getElementById('plots');
    var imagetype=document.getElementById('imagetype').value;
    var plot_str='';
    var category=document.getElementById('category');
    var metric=document.getElementById('metric');
    expansion_element=document.getElementById(arguments[0]);
    rows=document.getElementsByName(arguments[0]);
    if (expansion_element.innerHTML=='\u25B6'){
        expansion_element.innerHTML='\u25BC'
        show_row=arguments[0]+'_raw'+imagetype

        if (document.getElementById(show_row)==null){
            img=document.createElement('img')
            img.setAttribute('width',"600px")
            img.setAttribute('id',arguments[0]+'_raw'+imagetype)
            img.setAttribute('style','position:absolute;z-index:0')
            %s
            plots.appendChild(img)
        }else{
            document.getElementById(arguments[0]+'_raw'+imagetype).style.display=''
        }
        for (var i=0, il=rows.length;i<il;i++){
            rows[i].style.display='';
        }
    }else{
        expansion_element.innerHTML='\u25B6'
        document.getElementById(arguments[0]+'_raw'+imagetype).style.display='none'
        for (var i=0, il=rows.length;i<il;i++){
            rows[i].style.display='none';
        }
    }
}

function show_hide_categories(SelObject){
    var all_categories=document.getElementById('all_categories').value.split('$#!')
    var category=document.getElementById('category').value;
    var imagetype=document.getElementById('imagetype').value;
    var metric=document.getElementById('metric').value;
    for (var i=1, il=all_categories.length; i<il; i++){
        basename=metric+category+all_categories[i]
        raw_image=basename+'_raw'+imagetype
        ave_image=basename+'_ave'+imagetype
        checkbox=document.getElementsByName(raw_image)
        if (SelObject.value=='All'){
            if (checkbox[0].checked==false){
                checkbox[0].checked=true
                document.getElementById(ave_image).style.display=''
            }
        }else if (SelObject.value=='None'){
            if (checkbox[0].checked==true){
                checkbox[0].checked=false
                document.getElementById(ave_image).style.display='none'
            }
        }else if (SelObject.value=='Invert'){
            if (checkbox[0].checked==true){
                checkbox[0].checked=false
                document.getElementById(ave_image).style.display='none'
            }else if (checkbox[0].checked==false){
                checkbox[0].checked=true
                document.getElementById(ave_image).style.display=''
            }
        }
    }
    document.getElementById('show_category').selectedIndex=0;
}
</script>

</head>
<body>
<form action=''>
<input id="metric" type="hidden">
<input id="category" type="hidden">
<input id="imagetype" type="hidden" value="%s">
<input id="all_categories" type="hidden">
</form>
<table><tr>
<td><b>Select a Metric:</b></td>
<td>
<select onchange="javascript:changeMetric(this)" id="select_metric_combo">
<option>&nbsp;</option>
%s
</select>
</td>
<td><b>&nbsp;&nbsp;Select a Category:</b></td>
<td>
<select onchange="javascript:changeCategory(this)" id="select_category_combo">
<option>&nbsp;</option>
%s
</select>
</td>
</table>
<br>
<div style="width:950px">
<div id="plots" style="width:650px;height:550px;float:left;"></div>

<div id="legend" style="width:300px;height:550px;float:right;display:none;">
    <p><b>Show Categories:
    <select id="show_category" onchange="show_hide_categories(this);">
        <option value="">&nbsp;</option>
        <option value="All">All</option>
        <option value="None">None</option>
        <option value="Invert">Invert</option>
    </select>
    </b></p>
%s
<div style="position:relative;clear:both;">
<div style="position:relative;clear:both;display:none;" class="strong" id="nan_disclaimer">
<b>If the lines for some categories do not extend all the way to the right end of the x-axis, that means that at least one of the samples in that category does not have that many sequences.</b>
</div>
<br><br>
<table id="rare_data" border="1px">
%s
</table>
</div>
</div>
</body>
</html>
'''
