#!/usr/bin/env python
#file make_rarefaction_plots.py
from __future__ import division
__author__ = "Meg Pirrung"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Meg Pirrung", "Jesse Stombaugh"] 
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"

from matplotlib import use
use('Agg',warn=False)
from sys import exit
from qiime.parse import parse_rarefaction_data
from matplotlib.pylab import savefig,clf,gca,gcf,errorbar
import matplotlib.pyplot as plt
import os.path
from os.path import splitext, split
from qiime.colors import natsort, iter_color_groups
from qiime.util import create_dir
from numpy import isnan,nan,array,transpose,mean,std,arange

def save_ave_rarefaction_plots(xaxis, yvals, err, xmax, ymax, ops, \
            mapping_category, imagetype, res, data_colors, colors, fpath,\
            background_color,label_color,rarefaction_legend_mat,metric_name):
    '''This function creates the images, using matplotlib.'''
    
    #Create the plot image
    plt.clf()
    plt.title(metric_name + ": " + mapping_category)
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
    ax.set_ylabel("Rarefaction Measure: " + metric_name)
    
    #Create file for image
    imgpath = fpath+mapping_category+ '.'+imagetype
    
    #Save the image
    plt.savefig(imgpath, format=imagetype, dpi=res)

    #Get the image name for the saved image relative to the main directory
    image_loc = imgpath
    rarefaction_legend_mat[metric_name][mapping_category]['link']= \
                os.path.join('average_plots', metric_name + mapping_category + \
                             '.'+imagetype)
    
    return rarefaction_legend_mat

def save_single_rarefaction_plots(group,sample_dict,imagetype, metric_name,
                                  mapping_category, data_colors, colors,fpath,
                                  background_color,label_color,res,ymax,xmax,
                                  group_name,rarefaction_legend_mat):
    '''This function creates the images, using matplotlib.'''
    #Create the plot image
    plt.clf()
    plt.title(str(metric_name) + ": "+str(mapping_category)+": "+group_name)
    fig = plt.gcf()
    ax = fig.add_subplot(111)
    
    #Add the lines to the plot and create a dictionary to be used as the legend
    for o in group:
        if sample_dict.has_key(o):
            for i in sample_dict[o]:

                xaxis=[]
                #this creates duplicates of the xval, since there are several
                #iterations
                for t in range(len(sample_dict[o][i])):
                    xaxis.append(i)
                #Add values to the legend dictionary
                rarefaction_legend_mat[metric_name][mapping_category][group_name]['groupsamples'][o]=data_colors[colors[o]].toHex()
                
                #If all the yvals are nan at a particular xval, skip adding
                #it to the plot
                if not isnan(sample_dict[o][i])[0]:
                    scplot=ax.scatter(xaxis, sample_dict[o][i], \
                                        c=data_colors[colors[o]].toHex())

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
    ax.set_ylabel("Rarefaction Measure: " + str(metric_name))
     
    #Create file for image
    imgpath = fpath+ '.'+imagetype
    
    #Save the image
    plt.savefig(imgpath, format=imagetype, dpi=res)
    
    #Get the image name for the saved image relative to the main directory
    image_loc = imgpath
    
    #Add the image link to the legend dictionary

    rarefaction_legend_mat[metric_name][mapping_category][group_name]['link']=\
                os.path.join('all_other_plots',metric_name + mapping_category +\
                             group_name +  '.' + imagetype)

    return rarefaction_legend_mat

def get_rarefaction_data(rarefaction_data, col_headers):
    '''This function takes a rarefaction file and converts it into an array'''

    rare_mat_raw = array(rarefaction_data)
    rare_mat_min = [rare_mat_raw[x][2:] for x in range(0,len(rare_mat_raw))]
    seqs_per_samp = [rare_mat_raw[x][0] for x in range(0,len(rare_mat_raw))]
    sampleIDs = col_headers[3:]
    
    #Need to transpose the array to be used in averaging
    rare_mat_trans = transpose(array(rare_mat_min)).tolist()
    
    return rare_mat_trans, seqs_per_samp, sampleIDs

def ave_seqs_per_sample(matrix, seqs_per_samp, sampleIDs):
    """Calculate the average for each sampleID across each number of \
    seqs/sample"""
    
    ave_ser = {}
    temp_dict = {}
    #Iterate through the samples id's and create a dictionary
    for i,sid in enumerate(sampleIDs):
        temp_dict[sid] = {}
        for j,seq in enumerate(seqs_per_samp):
            try:
                temp_dict[sid][seq].append(matrix[i][j])
            except(KeyError):
                temp_dict[sid][seq] = []
                temp_dict[sid][seq].append(matrix[i][j])

    #create a dictionary for average data
    for sid in sampleIDs:
        ave_ser[sid] = []
        keys = temp_dict[sid].keys()
        keys.sort()
        for k in keys:
            ave_ser[sid].append(mean(array(temp_dict[sid][k]),0))
            
    return ave_ser

'''
#This function is currently not being used

def is_max_category_ops(mapping, mapping_category):
    """Count how many unique values there are for the supplied mapping \
    category and return true if all values are unique"""
    
    header = mapping[1]
    map_min = mapping[0]
    num_samples = len(map_min)
    index = header.index(mapping_category)
    seen = set()
    
    for m in map_min:
        seen.update([m[index]])
        
    return (len(seen) == num_samples), len(seen)
'''

def make_error_series(rare_mat, groups):
    """Create mean and error bar series for the supplied mapping category"""
    
    err_ser = dict()
    collapsed_ser = dict()
    seen = set()
    pre_err = {}
    
    ops = [k for k in groups]

    notfound = []             
    #Iterate through the groups
    for o in ops:
        pre_err[o] = []
        
        #For each sample in group, create a row in a list
        for samID in groups[o]:
            try:
                pre_err[o].append(rare_mat[samID])
            except(KeyError):
                notfound.append(samID)
                
        min_len = 1000 #1e10000
        #Iterate through the series data and convert it to float
        for series in pre_err[o]:
            series = [float(v) for v in series if v != 0]
            
            #determine the minimum length of a series
            if len(series) < min_len:
                min_len = len(series)
        
        pre_err[o] = [x[:min_len] for x in pre_err[o]]
    
    #iterate through the groups and calculate std deviations and error
    for o in ops:
        opsarray = array(pre_err[o])
        mn = mean(opsarray, 0)
        collapsed_ser[o] = mn.tolist()
        stddev = std(opsarray, 0)
        err_ser[o] = stddev.tolist()
        
    return collapsed_ser, err_ser, ops

'''
function is not used

def get_overall_averages(rare_mat, sampleIDs):
    """Make series of averages of all values of seqs/sample for each \
    sampleID"""
    
    overall_ave = dict();
    for s in sampleIDs:
        overall_ave[s] = mean(array(rare_mat[s]))
        
    return overall_ave
'''

def save_rarefaction_data(rare_mat,xaxis, xmax, \
mapping_category, colors, rare_type, data_colors, groups):
    '''This function formats the average data and writes it to the output 
       directory'''
     
    #get the error data
    yaxis, err, ops = make_error_series(rare_mat,groups)
   
    lines = []
    lines.append("# "+rare_type+'\n')
    lines.append("# "+mapping_category+'\n')
    line = ''
    line += 'xaxis: '
    for v in xaxis:
        line += str(v) + '\t'
    line += '\n'
    lines.append(line)
    lines.append('xmax: '+str(xmax)+'\n')
    
    for o in colors.keys():
        lines.append(">> " + o + '\n')
        
        #write the color lines
        if colors != None:
            try:
                lines.append("color " + data_colors[colors[o]].toHex() + '\n')
            except(KeyError):
                print 'Color reference is missing!'
        
        #write the rarefection series lines
        lines.append('series ')
        line = ''
        try:
            for v in yaxis[o]:
                line += str(v) + '\t'
        except(TypeError):
            line += str(yaxis[o])
        line += '\n'
        lines.append(line)
        
        #write the rarefaction error lines
        lines.append('error ')
        line = ''
        try:
            for e in err[o]:
                line += str(e) + '\t'
        except(TypeError):
            line += str(err[o])
        line += '\n'
        lines.append(line)
    
    return lines   

def make_averages(color_prefs, data, background_color, label_color, rares, \
                    output_dir,resolution,imagetype):
    '''This is the main function, which takes the rarefaction files, calls the
        functions to make plots and formatting the output html.''' 
    rarelines = []
    rarefaction_legend_mat={}
    
    #Create the directories, where plots and data will be written
    all_output_dir = os.path.join(output_dir, 'all_other_plots')
    create_dir(all_output_dir)
    ave_output_dir = os.path.join(output_dir, 'average_plots')
    create_dir(ave_output_dir)
    ave_data_file_path=os.path.join(output_dir,'average_tables')
    create_dir(ave_data_file_path,False)
    
    metric_num=0
    rarefaction_legend_mat={}
    rarefaction_data_mat={}
    rare_num=0
    
    #Iterate through the rarefaction files
    for r in natsort(rares):
        raredata = rares[r]
        metric_name=r.split('.')[0]
        
        #convert the rarefaction data into variables
        col_headers,comments,rarefaction_fn,rarefaction_data=rares[r]
        
        #Here we only need to perform these steps once, since the data is
        #the same for all rarefaction files
        if rare_num==0:
            
            #Remove samples from the mapping file, which contain no data after
            #rarefaction
            updated_mapping=[]
            for j in data['map']:
                
                #Add the mapping header
                if j[0]=='SampleID':
                    updated_mapping.append(j)
                
                #Determine if the sample exists in the rarefaction file
                for i in col_headers[3:]:
                    if j[0]==i:
                        updated_mapping.append(j)
            
            #Get the groups and colors for the updated mapping file
            groups_and_colors=iter_color_groups(updated_mapping,color_prefs)
            groups_and_colors=list(groups_and_colors)
        
        #parse the rarefaction data
        
        rare_mat_trans, seqs_per_samp, sampleIDs = \
        get_rarefaction_data(rarefaction_data, col_headers)
        
        rarefaction_legend_mat[metric_name]={}
        
        #Create dictionary variables and get the colors for each Sample
        sample_colors=None
        for i in range(len(groups_and_colors)):
            labelname=groups_and_colors[i][0]
            
            #Create a legend dictionary for html output
            rarefaction_legend_mat[metric_name][labelname]={}
            
            #If this is the first time iterating through the rarefaction data
            #create a data dictionary for html output
            if rare_num==0:
                rarefaction_data_mat[labelname]={}
   
            #If the labelname is SampleID, use the colors assigned
            if labelname=='SampleID':
                sample_colors=groups_and_colors[i][2]
                sample_data_colors=groups_and_colors[i][3]
        
        rare_num=1
        
        #If sample colors were not assigned, create a list of sample colors
        if not sample_colors:
            samples_and_colors=iter_color_groups(updated_mapping, \
                {'SampleID': {'column': 'SampleID', 'colors': \
                (('red', (0, 100, 100)), ('blue', (240, 100, 100)))}})
            samples_and_colors=list(samples_and_colors)
            sample_colors=samples_and_colors[0][2]
            sample_data_colors=samples_and_colors[0][3]
        
        sample_dict = {}
        
        #Create a dictionary containing the samples
        for i,sid in enumerate(sampleIDs):
            sample_dict[sid] = {}
            for j,seq in enumerate(seqs_per_samp):
                try:
                    sample_dict[sid][seq].append(rare_mat_trans[i][j])
                except(KeyError):
                    sample_dict[sid][seq] = []
                    sample_dict[sid][seq].append(rare_mat_trans[i][j])

        #convert xvals to float
        xaxisvals = [float(x) for x in set(seqs_per_samp)]
        xaxisvals.sort()
        
        #get the rarefaction averages
        rare_mat_ave = ave_seqs_per_sample(rare_mat_trans, seqs_per_samp, \
        sampleIDs)

        #calculate the max xval
        xmax = max(xaxisvals) + (xaxisvals[len(xaxisvals)-1] - \
        xaxisvals[len(xaxisvals)-2])
        
        '''
        #get the overall average
        #overall_average = get_overall_averages(rare_mat_ave, sampleIDs)
        

        rarelines.append("#" + r + '\n')
          
        for s in sampleIDs:
            rarelines.append('%f'%overall_average[s] + '\n')
        '''
        
        #iterate through the groups
        for i in range(len(groups_and_colors)):
            labelname=groups_and_colors[i][0]
            groups=groups_and_colors[i][1]
            colors=groups_and_colors[i][2]
            data_colors=groups_and_colors[i][3]
            data_color_order=groups_and_colors[i][4]
            ave_file_path=os.path.join(ave_data_file_path,metric_name)
            
            #save the rarefaction averages
            
            rare_lines=save_rarefaction_data(rare_mat_ave, xaxisvals, xmax, \
                                    labelname, colors, r, data_colors, groups)
            
            #write out the rarefaction average data
            open(ave_file_path+labelname+'.txt','w').writelines(rare_lines)
            
            #take the formatted rarefaction averages and format the results
            rares_data = parse_rarefaction_data( \
                                        ''.join(rare_lines[:]).split('\n'))
            
            #determine the ymax based on the average data
            #multiple the ymax, since the dots can end up on the border
            ymax = (max([max(v) for v in rares_data['series'].values()]) + \
                    max([max(v) for v in rares_data['error'].values()]))*1.15
            
            #Iterate through the groups and create the legend dictionary
            for g in groups:
                #generate the filepath for the image file
                file_path = os.path.join(all_output_dir, \
                                            metric_name+labelname+g)
                #create a dictionary of samples and their colors
                rarefaction_legend_mat[metric_name][labelname][g]={}
                rarefaction_legend_mat[metric_name][labelname][g]['groupsamples']={}
                rarefaction_legend_mat[metric_name][labelname][g]['groupcolor']=\
                                        data_colors[colors[g]].toHex()
                
                #Create raw plots for each group in a category
                rarefaction_legend_mat=save_single_rarefaction_plots( \
                                            groups[g],sample_dict, \
                                            imagetype,metric_name,labelname, \
                                            sample_data_colors,sample_colors, \
                                            file_path,background_color, \
                                            label_color,resolution,ymax,xmax,g, \
                                            rarefaction_legend_mat)
            
            #Create the average rarefaction plots
            rarefaction_data_mat,rarefaction_legend_mat=make_plots(\
                                background_color, label_color, \
                                rares_data, ymax, xmax,ave_output_dir, \
                                resolution, imagetype,groups, colors, \
                                data_colors,metric_name,labelname, \
                                rarefaction_data_mat,rarefaction_legend_mat)

    #format the html output
    html_output=make_html(rarefaction_legend_mat, \
                            rarefaction_data_mat,xaxisvals,imagetype)
    
    return html_output
    
def make_html(rarefaction_legend_mat,rarefaction_data_mat,xaxisvals,imagetype):
    
    legend_td=[]
    summarized_table=[]
    metric_select_html=[]
    category_select_html=[]
    data_table_html=[]
    metrics=[]
    category_colors={}
    cat_iter=0
    
    #iterate the legend dictionary
    for m in natsort(rarefaction_legend_mat):
        #Create the metric select box options
        metric_select_html.append('<option value="%s">%s</option>' % (m,m))
        metrics.append(m)
        
        #iterate through the categories in the legend dictionary
        for category in natsort(rarefaction_legend_mat[m]):
            if cat_iter==0:
                #Create the select box options
                category_select_html.append('<option value="%s">%s</option>' % \
                                (category,category))
            
            #iterate through the groups in the legend dictionary and create
            #the html formatted rows for each category and group
            for group in natsort(rarefaction_legend_mat[m][category]):
                if group <> 'link':
                   
                    category_colors[group]=\
                       rarefaction_legend_mat[m][category][group]['groupcolor']
                    legend_td.append('<tr name="%s" style="display: none;" onmouseover="javascript:change_plot(\'./%s\')" onmouseout="javascript:change_plot(\'./%s\')"><td class="data" bgcolor="%s"><b>%s</b></td></tr>' % (m+category,rarefaction_legend_mat[m][category][group]['link'],rarefaction_legend_mat[m][category]['link'],rarefaction_legend_mat[m][category][group]['groupcolor'], group)) 

        cat_iter=1
        
    #iterate through the data dictionary and format the rows for the html 
    #data table
    for category in rarefaction_data_mat:
        data_table_html.append('<tr name="%s" style="display: none;"><td class="headers">%s</td><td class="headers">Seqs/Sample</td>' % (category,category))
        for j in metrics:
            data_table_html.append('<td class="headers">%s Ave.</td><td class="headers">%s Err.</td>' % (j,j))
        data_table_html.append('</tr>')
        data_table_html.append('<tr name="%s" style="display: none;"></tr>' % (category))
        for g in natsort(rarefaction_data_mat[category]):
            for i in range(len(xaxisvals)):
                
                data_table_html.append('<tr name="%s" style="display: none;">' % (category))
                data_table_html.append('<td class="data" bgcolor="%s">%s</td><td class="data">%s</td>' % (category_colors[g],g,xaxisvals[i]))
                for m in rarefaction_data_mat[category][g]:
                    data_table_html.append('<td class="data">%s</td><td class="data">%s</td>' % (rarefaction_data_mat[category][g][m]['ave'][i],rarefaction_data_mat[category][g][m]['err'][i]))
        data_table_html.append('</tr>')
    
    #Create the table that contains the plots and table  
    plot_html='<table><tr><td><img id="plots" style="display: none; height: 400px;" \></td><td id="legend" style="display: none;"><b>Legend</b><div STYLE="border: thin black solid; height: 300px; width: 150px; font-size: 12px; overflow: auto;"><table>%s</table></div></td></tr></table>' % \
                            ('\n'.join(legend_td))
    
    #insert the formatted rows into the html string at the bottom of this file 
    html_output=HTML % ('.'+imagetype,
                        '\n'.join(metric_select_html), \
                        '\n'.join(category_select_html), \
                        plot_html, \
                        '\n'.join(data_table_html))

    return html_output
    
def make_plots(background_color, label_color, rares, ymax, xmax,\
                output_dir, resolution, imagetype,groups,colors,data_colors, \
                metric_name,labelname,rarefaction_data_mat,\
                rarefaction_legend_mat):   
    '''This is the main function for generating the rarefaction plots and html
        file.'''
    
    #set the arrays for html text.
    plot_table_html=[]
    data_table_html=[]
    category_select_html=[]
    legend_html=[]
    legend_td=[]
    
    #Get the alpha rare data
    raredata = rares
    
    #generate the filepath for the image file
    file_path = os.path.join(output_dir, \
    splitext(split(raredata['headers'][0])[1])[0])

    #Sort and iterate through the groups
    for i in natsort(groups):

        #for k in groups[i]:
        for j in range(len(raredata['xaxis'])):
            group_field=i
            
            seq_per_sample_field=int(raredata['xaxis'][j])
            color_field=data_colors[colors[group_field]].toHex()
             
            #If a field is missing, then it means that one of the 
            #samples did not contain enough sequences.
            #For this case, we will assign the value as n.a.
            try:
                average_field=raredata['series'][i][j]
                error_field=raredata['error'][i][j]
                if isnan(average_field):
                    error_field=nan
            except:
                average_field=nan
                error_field=nan
            
            #Add context to the data dictionary, which will be used in the html
            if rarefaction_data_mat[labelname].has_key(i):
                if rarefaction_data_mat[labelname][i].has_key(metric_name):
                    rarefaction_data_mat[labelname][i][metric_name]['ave'].append(''.join('%10.3f' % (raredata['series'][i][j])))
                    rarefaction_data_mat[labelname][i][metric_name]['err'].append(''.join('%10.3f' % (raredata['error'][i][j])))
                else:
                    rarefaction_data_mat[labelname][i][metric_name]={}
                    rarefaction_data_mat[labelname][i][metric_name]['ave']=[]
                    rarefaction_data_mat[labelname][i][metric_name]['err']=[]
                    rarefaction_data_mat[labelname][i][metric_name]['ave'].append(''.join('%10.3f' % (raredata['series'][i][j])))
                    rarefaction_data_mat[labelname][i][metric_name]['err'].append(''.join('%10.3f' % (raredata['error'][i][j])))
            else:
                rarefaction_data_mat[labelname][i]={}
                rarefaction_data_mat[labelname][i][metric_name]={}
                rarefaction_data_mat[labelname][i][metric_name]['ave']=[]
                rarefaction_data_mat[labelname][i][metric_name]['err']=[]
                rarefaction_data_mat[labelname][i][metric_name]['ave'].append(''.join('%10.3f' % (raredata['series'][i][j])))
                rarefaction_data_mat[labelname][i][metric_name]['err'].append(''.join('%10.3f' % (raredata['error'][i][j])))
            
    categories = [k for k in groups]

    #Create the rarefaction average plot and get updated legend information
    rarefaction_legend_mat=\
    save_ave_rarefaction_plots(raredata['xaxis'], raredata['series'], \
                           raredata['error'], xmax, ymax, categories, \
                           labelname, imagetype, resolution, data_colors, \
                           colors, file_path, background_color, label_color, \
                           rarefaction_legend_mat, metric_name)
    
    return rarefaction_data_mat,rarefaction_legend_mat
    
    
HTML='''
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
  <meta http-equiv="content-type" content="text/html;">
  <title>Rarefaction Curves</title>
<style type="text/css">
td.data{font-size:10px}
td.headers{font-size:12px}
table{border-spacing:0px}
</style>
<script language="javascript" type="text/javascript">

function change_plot(object){
plot_img=document.getElementById('plots');
plot_img.src=object
}


function changeMetric(SelObject){
    var old_category=document.getElementById('category');
    var old_metric=document.getElementById('metric');
    var imagetype=document.getElementById('imagetype').value;
    var legend=document.getElementById('legend');

    if (old_category.value != ''){
        plot_img=document.getElementById('plots');
        plot_loc='./average_plots/'+SelObject.value+old_category.value+imagetype
        plot_img.src=plot_loc
        plot_img.style.display="";
        legend.style.display="";
         data_display=document.getElementsByName(SelObject.value+old_category.value)
        for (var i=0;i<data_display.length;i++){
            data_display[i].style.display="";
        }
        data_hide=document.getElementsByName(old_metric.value+old_category.value)
        for (var i=0;i<data_hide.length;i++){
            data_hide[i].style.display="none";
        }
        
        data_display=document.getElementsByName(old_category.value)
        for (var i=0;i<data_display.length;i++){
            data_display[i].style.display="";
        }
    }
    
old_metric.value=SelObject.value;
}

function changeCategory(SelObject){
    var old_category=document.getElementById('category');
    var old_metric=document.getElementById('metric');
    var imagetype=document.getElementById('imagetype').value;
    var legend=document.getElementById('legend');

    if (old_metric.value != ''){
        plot_img=document.getElementById('plots');
        plot_loc='./average_plots/'+old_metric.value+SelObject.value+imagetype;
        plot_img.src=plot_loc;
        plot_img.style.display="";
        legend.style.display="";
        
        data_display=document.getElementsByName(SelObject.value)
        for (var i=0;i<data_display.length;i++){
            data_display[i].style.display="";
        }

        data_hide=document.getElementsByName(old_category.value)
        for (var i=0;i<data_hide.length;i++){
            data_hide[i].style.display="none";
        }

        data_display=document.getElementsByName(old_metric.value+SelObject.value)
        for (var i=0;i<data_display.length;i++){
            data_display[i].style.display="";
        }

        data_hide=document.getElementsByName(old_metric.value+old_category.value)
        for (var i=0;i<data_hide.length;i++){
            data_hide[i].style.display="none";
        }
    }
old_category.value=SelObject.value;
}

</script>
</head>
<body>
<form>
<input id="metric" type="hidden"></input>
<input id="category" type="hidden"></input>
<input id="imagetype" type="hidden" value="%s"></input>
</form>
<table><tr>
<td><b>Select a Metric:</b></td>
<td>
<select onchange="javascript:changeMetric(this)">
<option></option>
%s
</select>
</td>
<td><b>Select a Category:</b></td>
<td>
<select onchange="javascript:changeCategory(this)">
<option></option>
%s
</select>
</td>
</table>
<br>
%s
<br>
<table id="rare_data">
%s
</table>
</body>
</html>
'''
    