#!/usr/bin/env python
# File created Sept 29, 2010
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Release"
 
from cogent.parse.fasta import MinimalFastaParser
from numpy import arange, std, average
from pylab import plot, savefig, xlabel, ylabel, text, \
    hist, figure, legend, title, show, xlim ,ylim, xticks, yticks,\
    scatter, subplot
from matplotlib.font_manager import fontManager, FontProperties

from qiime.parse import parse_qual_score



    
        
        
def bin_qual_scores(qual_scores):
    """ Bins qual score according to nucleotide position 
    
    qual_scores: Dict of label: numpy array of base scores
    """
    
    qual_bins = []
    
    qual_lens = []
    
    for l in qual_scores.values():
        qual_lens.append(len(l))
        
    max_seq_size = max(qual_lens)
    
    for base_position in range(max_seq_size):
        qual_bins.append([])
        
        for scores in qual_scores.values():
            # Add score if exists in base position, otherwise skip
            try:
                qual_bins[base_position].append(scores[base_position])
            except IndexError:
                continue
    
    
    return qual_bins
    
def get_qual_stats(qual_bins):
    """ Generates bins of averages, std devs, total NT from quality bins"""
    
    ave_bins = []
    std_dev_bins = []
    total_bases_bins = []
    
    for base_position in qual_bins:
        
        total_bases_bins.append(len(base_position))
        
        std_dev_bins.append(std(base_position))
        
        ave_bins.append(average(base_position))
    
    
    return ave_bins, std_dev_bins, total_bases_bins
    
    
    
def plot_qual_report(ave_bins,
                     std_dev_bins,
                     total_bases_bins,
                     score_min,
                     output_dir):
    """ Plots, saves graph showing quality score averages, stddev.
    
    Additionally, the total nucleotide count for each position is shown on 
     a second subplot
     
    ave_bins: list with average quality score for each base position
    std_dev_bins: list with standard deviation for each base position
    total_bases_bins: list with total counts of bases for each position
    score_min: lowest value that a given base call can be and still be
     acceptable.  Used to generate a dotted line on the graph for easy assay
     of the poor scoring positions.
    output_dir: output directory
    """
    
    t = arange(0, len(ave_bins), 1)
    
    std_dev_plus = []
    std_dev_minus = []
    
    for n in range(len(ave_bins)):
        std_dev_plus.append(ave_bins[n] + std_dev_bins[n])
        std_dev_minus.append(ave_bins[n] - std_dev_bins[n])
    
    figure_num = 0
    
    f = figure(figure_num, figsize=(8, 10))
    
    figure_title = "Quality Scores Report"
    
    f.text(.5, .93, figure_title, horizontalalignment='center', size="large")
    
    subplot(2, 1, 1)
    
    plot(t, ave_bins, linewidth=2.0, color="black")
    plot(t, std_dev_plus, linewidth=0.5, color="red")
    
    dashed_line = [score_min] * len(ave_bins)
    
    l, = plot(dashed_line, '--', color='gray')
    
    plot(t, std_dev_minus, linewidth=0.5, color="red")
    
    

    
    legend( ('Quality Score Average', 'Std Dev', 'Score Threshold') , loc = 'lower left')

    
    xlabel("Nucleotide Position")
    ylabel("Quality Score")
    
    
    subplot(2, 1, 2)
    
    plot(t, total_bases_bins, linewidth=2.0, color="blue")
    
    xlabel("Nucleotide Position")
    ylabel("Nucleotide Counts")
    
    
    outfile_name = output_dir + "/quality_scores_plot.pdf"
    
    savefig(outfile_name)

    
def write_qual_report(ave_bins,
                      std_dev_bins,
                      total_bases_bins,
                      output_dir):
    """ Writes data in bins to output text file
    
    ave_bins: list with average quality score for each base position
    std_dev_bins: list with standard deviation for each base position
    total_bases_bins: list with total counts of bases for each position
    output_dir: output directory
    """
    
    outfile_name = output_dir + "/quality_bins.txt"
    
    outfile = open(outfile_name, "w")
    
    outfile.write("# Average quality score bins\n")
    
    outfile.write(",".join(str("%2.3f" % ave) for ave in ave_bins) + "\n")
    
    outfile.write("# Standard deviation bins\n")
    
    outfile.write(",".join(str("%2.3f" % std) for std in std_dev_bins) + "\n")
    
    outfile.write("# Total bases per nucleotide position bins\n")
    
    outfile.write(",".join(str("%d" %\
     total_bases) for total_bases in total_bases_bins))
    
def generate_histogram(qual_fp,
                       output_dir,
                       score_min=25):
    """ Main program function for generating quality score histogram

    qual_fp: quality score filepath
    output_dir: output directory
    score_min: minimum score to be considered a reliable base call, used 
     to generate dotted line on histogram for easy visualization of poor
     quality scores."""
    
    qual_lines = open(qual_fp, "U")
    
    qual_scores = parse_qual_score(qual_lines)
    
   
    # Sort bins according to base position
    qual_bins = bin_qual_scores(qual_scores)
    
    # Get average, std dev, and total nucleotide counts for each base position
    ave_bins, std_dev_bins, total_bases_bins = get_qual_stats(qual_bins)
    
    plot_qual_report(ave_bins, std_dev_bins, total_bases_bins, score_min,
     output_dir)
     
    # Save values to output text file
    write_qual_report(ave_bins, std_dev_bins, total_bases_bins, output_dir)
