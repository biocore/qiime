# Runs vegan function Moran.I on QIIME distance matrix
# usage:
# R --slave --args --source_dir $QIIME_HOME/qiime/support_files/R/ -d unifrac.txt -m Fasting_Map.txt -c Treatment -o morans_i < morans_i.r
#
# print help string:
# R --slave --args -h --source_dir $QIIME_HOME/qiime/support_files/R/ < morans_i.r
#
# Requires command-line param --source_dir pointing to QIIME R source dir

# load libraries and source files
args <- commandArgs(trailingOnly=TRUE)
if(!is.element('--source_dir', args)){
    stop("\n\nPlease use '--source_dir' to specify the R source code directory.\n\n")
}
sourcedir <- args[which(args == '--source_dir') + 1]
source(sprintf('%s/loaddata.r',sourcedir))
source(sprintf('%s/util.r',sourcedir))
load.library('optparse')
load.library('ape')

# make option list and parse command line
option_list <- list(
    make_option(c("--source_dir"), type="character",
        help="Path to R source directory [required]."),
    make_option(c("-d", "--distmat"), type="character",
        help="Input distance matrix [required]."),
    make_option(c("-m", "--mapfile"), type="character",
        help="Input metadata mapping file [required]."),
    make_option(c("-c", "--category"), type="character",
        help="Metadata column header giving cluster IDs [required]"),
    make_option(c("-o", "--outdir"), type="character", default='.',
        help="Output directory [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# error checking
if(is.null(opts$mapfile)) stop('Please supply a mapping file.')
if(is.null(opts$category)) stop('Please supply a mapping file header.')
if(is.null(opts$distmat)) stop('Please supply a distance matrix.')

# create output directory if needed
if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)

# load qiime data
map <- load.qiime.mapping.file(opts$mapfile)
distmat <- load.qiime.distance.matrix(opts$distmat)
qiime.data <- remove.nonoverlapping.samples(map=map, distmat=distmat)

# run s
# use inverse distance as weight. Set all zeros to very small number, otherwise
# Moran.I will error out.
qiime.data$distmat[qiime.data$distmat == 0] <- .Machine$double.eps
weights = 1/qiime.data$distmat

# set diagonal back to zero
diag(weights) <- 0

results <- Moran.I(x=qiime.data$map[[opts$category]],weight=weights)

# write output file
filepath <- sprintf('%s/morans_i_results.txt',opts$outdir)
sink(filepath)
print(results)
sink(NULL)
