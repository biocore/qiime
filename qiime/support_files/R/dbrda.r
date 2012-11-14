# Runs vegan function capscale (db-RDA) on QIIME distance matrix
# usage:
# R --slave --args --source_dir $QIIME_HOME/qiime/support_files/R/ -d
# unifrac.txt -m Fasting_Map.txt -c Treatment -o dbrda < dbrda.r
#
# print help string:
# R --slave --args -h --source_dir $QIIME_HOME/qiime/support_files/R/ < dbrda.r
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
load.library('vegan')

# Make option list and parse command line.
option_list <- list(
    make_option(c("--source_dir"), type="character",
        help="Path to R source directory [required]."),
    make_option(c("-d", "--distmat"), type="character",
        help="Input distance matrix [required]."),
    make_option(c("-m", "--mapfile"), type="character",
        help="Input metadata mapping file [required]."),
    make_option(c("-c", "--category"), type="character",
        help="Metadata column header giving cluster IDs [required]"),
    make_option(c("-n", "--num_permutations"), type="integer", default=999,
        help="Number of permutations [default %default]."),
    make_option(c("-o", "--outdir"), type="character", default='.',
        help="Output directory [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Make sure we have our required files.
if (is.null(opts$mapfile)) stop('Please supply a mapping file.')
if (is.null(opts$category)) stop('Please supply a mapping file category.')
if (is.null(opts$distmat)) stop('Please supply a distance matrix.')

# Create output directory if needed.
if (opts$outdir != ".") dir.create(opts$outdir, showWarnings=FALSE, recursive=TRUE)

# Load data.
map <- load.qiime.mapping.file(opts$mapfile)
distmat <- load.qiime.distance.matrix(opts$distmat)
qiime.data <- remove.nonoverlapping.samples(map=map, distmat=distmat)

# More error checking.
if (nrow(qiime.data$map) == 0)
    stop('\n\nMapping file and distance matrix have no samples in common.\n\n')
if (!is.element(opts$category, colnames(qiime.data$map)))
    stop(sprintf('\n\nHeader %s not found in mapping file.\n\n', opts$category))

# Run db-RDA and create a plot of the results.
factor = as.factor((qiime.data$map[[opts$category]]))
factors.frame <- data.frame(factor)
capscale.results <- capscale(as.dist(qiime.data$distmat) ~ factor, factors.frame)

# Perform significance test.
sig.test.results <- permutest(capscale.results,
                              permutations=opts$num_permutations)

# Write out results files.
capscale.results.filepath <- sprintf('%s/dbrda_results.txt', opts$outdir)
sink(capscale.results.filepath)
print(capscale.results)
print(sig.test.results)
sink(NULL)

plot.filepath <- sprintf('%s/dbrda_plot.pdf', opts$outdir)
pdf(plot.filepath)
plot(capscale.results, display=c("wa", "bp"))
