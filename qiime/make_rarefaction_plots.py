#!/usr/bin/env python
#file make_rarefaction_plots.py
from __future__ import division
__author__ = "Meg Pirrung"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Meg Pirrung", "Jesse Stombaugh"] 
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Meg Pirrung"
__email__ = "meg.pirrung@colorado.edu"
__status__ = "Development"

"""
Author: Meg Pirrung (meg.pirrung@colorado.edu) 
Status: Prototype

Requirements:
Python 2.5
Matplotlib
Numpy
"""
from matplotlib import use,rc
use('Agg',warn=False)
from string import strip
from matplotlib.pylab import savefig, clf, gca, gcf, errorbar
import matplotlib.pyplot as plt
import os.path
from os.path import exists, splitext, split
from qiime.colors import Color, natsort, iter_color_groups
from qiime.util import create_dir
from numpy import isnan
def save_rarefaction_plots(xaxis, yvals, err, xmax, ymax, ops, \
            mapping_category, itype, res, rtype, data_colors, colors, fpath,\
            background_color,label_color):
    '''This function creates the images, using matplotlib.'''
    
    #Creathe plot image
    plt.clf()
    plt.title((splitext(split(rtype)[1])[0]) + ": " + mapping_category)
    fig  = plt.gcf()
     
    #Add the lines to the plot
    for o in ops:
        l = o
        plt.errorbar(xaxis[:len(yvals[o])], yvals[o], \
                     yerr=err[o][:len(yvals[o])], label = l, color = \
                     data_colors[colors[o]].toHex(),elinewidth=1,lw=2,capsize=4)
    
    #get the plot axis
    ax = plt.gca()
    ax.set_axis_bgcolor(background_color)
    
    #set tick colors and width
    for line in ax.yaxis.get_ticklines():
        # line is a matplotlib.lines.Line2D instance
        line.set_color(label_color)
        line.set_markeredgewidth(1)
        
    for line in ax.xaxis.get_ticklines():
        # line is a matplotlib.lines.Line2D instance
        line.set_color(label_color)
        line.set_markeredgewidth(1)
    
    #set x/y limits and labels for plot
    ax.set_axisbelow(True)
    ax.set_xlim((0,xmax))
    ax.set_ylim((0,ymax))
    ax.set_xlabel('Sequences Per Sample')
    ax.set_ylabel("Rarefaction Measure: " + splitext(split(rtype)[1])[0])
    
    #Create file for image
    imgpath = os.path.join(fpath,mapping_category) + '.'+itype
    
    #Save the image
    plt.savefig(imgpath, format=itype, dpi=res)

    #Get the image name for the saved image relative to the main directory
    image_loc = splitext(split(rtype)[1])[0] + "/" + mapping_category +"." + \
                          itype

    return image_loc
    
def make_plots(color_prefs, data, background_color, label_color, rares, ymax, \
                output_dir, resolution, imagetype):   
    '''This is the main function for generating the rarefaction plots and html
        file.'''
    
    #set the arrays for html text.
    plot_table_html=[]
    data_table_html=[]
    category_select_html=[]
    legend_html=[]
    
    #Assign colors based on the supplied groups
    groups_and_colors=iter_color_groups(data['map'],color_prefs)
    groups_and_colors=list(groups_and_colors)

    #Iterate through the collated alpha files.
    iterator=0
    for r in natsort(rares):
        
        #Get the Category Name
        group_name=r.split('.')[0]
        
        #Check that the group is present in the prefs file.
        if color_prefs.has_key(group_name):
            
            #Get the alpha rare data
            raredata = rares[r]
            legend_td=['<td class="headers"><b>%s Coloring:</b></td>' % (group_name)]
            
            #generate the filepath for the image file
            file_path = os.path.join(output_dir, \
            splitext(split(raredata['headers'][0])[1])[0])
        
            #If the folder is not present, create it.
            if not os.path.exists(file_path):
                create_dir(file_path,False)
            if ymax == 0:
                ymax = max([max(v) for v in raredata['series'].values()]) + \
                max([max(v) for v in raredata['error'].values()])
            xmax = raredata['xaxis'][len(raredata['xaxis'])-1]
            
            groups=groups_and_colors[iterator][1]
            colors=groups_and_colors[iterator][2]
            data_colors=groups_and_colors[iterator][3]

            #Sort and iterate through the groups
            for i in natsort(groups):
                #Create the legend rows
                legend_td.append('<td class="data" bgcolor="%s">%s</td>' % (data_colors[colors[i]].toHex(), i))

                
                for k in groups[i]:
                    for j in range(len(raredata['xaxis'])):
                        group_field=i
                        sample_field=k
                        seq_per_sample_field=int(raredata['xaxis'][j])
                        color_field=data_colors[colors[i]].toHex()
                        
                        #If a field is missing, then it means that one of the 
                        #samples did not contain enough sequences.
                        #For this case, we will assign the value as n.a.
                        try:
                            average_field=raredata['series'][i][j]
                            error_field=raredata['error'][i][j]
                            if isnan(average_field):
                                error_field='nan'
                        except:
                            average_field='nan'
                            error_field='nan'
                            
                        #Create the data table rows
                        data_table_html.append('<tr name="%s" style="display: none;"><td class="data" bgcolor="%s">%s</td><td class="data">%s</td><td class="data">%s</td><td class="data">%s</td><td class="data">%s</td></tr>' \
                            % (group_name, color_field, group_field, \
                               sample_field, seq_per_sample_field, \
                               average_field, error_field))
                    
            #Add to the legend rows.
            legend_html.append('<tr name="%s" style="display: none;">%s</tr>' \
                                % (group_name,'\n'.join(legend_td)))
            
            categories = [k for k in groups]
            
            #Create the image and return the file location
            image_loc=\
            save_rarefaction_plots(raredata['xaxis'], raredata['series'], \
                                   raredata['error'], xmax, ymax, categories, \
                                   raredata['headers'][1],  imagetype, \
                                   resolution, raredata['headers'][0], \
                                   data_colors, colors, file_path, \
                                   background_color, label_color)
            
            #Write the image plot rows
            plot_table_html.append('<tr name="%s" style="display: none;"><td colspan="20"><img src="./%s" \></td></tr>' % \
                                    (group_name, image_loc))
            
            #Create the select box options
            category_select_html.append('<option value="%s">%s</option>' % \
                                        (group_name,group_name))
                                    
                                    
            iterator=iterator+1
    
    
        #Generate the html output
        html_output=HTML % ('\n'.join(category_select_html), \
                              '\n'.join(plot_table_html), \
                              '\n'.join(legend_html), \
                              '\n'.join(data_table_html))

    return html_output

    
HTML='''
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
  <meta http-equiv="content-type" content="text/html;">
  <title>Rarefaction Curves</title>
<style type="text/css">
td.data{font-size:10px}
td.headers{font-size:12px}
</style>
<script language="javascript" type="text/javascript">

function changeCategory(SelObject){

var header_display=document.getElementById('rare_header');
header_display.style.display='';

var old_selected=document.getElementById('old_selected');

data_display=document.getElementsByName(SelObject.value)
for (var i=0;i<data_display.length;i++){
    data_display[i].style.display="";
}

data_hide=document.getElementsByName(old_selected.value)
for (var i=0;i<data_hide.length;i++){
    data_hide[i].style.display="none";
}

old_selected.value=SelObject.value;
}
</script>
</head>
<body>
<form>
<input id="old_selected" type="hidden" value="" \>
</form>
<table><tr>
<td><b>Select a Category:</b></td>
<td>
<select onchange="javascript:changeCategory(this)">
<option></option>
%s
</select>
</td>
</table>
<br>
<table id="rare_plots">
%s
</table>
<br>
<table id="legend">
%s
</table>
<br>
<table id="rare_data">
<tr id="rare_header" style="display: none;">
<td class="headers"><b>Category</b></td>
<td class="headers"><b>SampleID</b></td>
<td class="headers"><b>Sequence Depth</b></td>
<td class="headers"><b>Rarefaction Ave.</b></td>
<td class="headers"><b>Error</b></td>
</tr>
%s
</table>
</body>
</html>
'''
    