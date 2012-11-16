# Runs detrending on QIIME pcoa table
#
# standalone usage:
# run with R --vanilla --slave --args -p pcfile -m mapfile -c category -o outdir --source_dir $QIIME_HOME/qiime/support_files/R/ < detrend.r
# 
# print help string:
# R --slave --args -h --source_dir $QIIME_HOME/qiime/support_files/R/ < detrend.r
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
source(sprintf('%s/detrend_util.r',sourcedir))
load.library('optparse')
load.library('RColorBrewer')

# make option list and parse command line
option_list <- list(
    make_option(c("--source_dir"), type="character",
        help="Path to R source directory [required]."),
    make_option(c("-i","--pcoa"), type="character",
        help="Path to PCOA file [required]."),
    make_option(c("-m", "--mapfile"), type="character",
        help="Input metadata mapping file [optional]."),
    make_option(c("-c", "--category"), type="character",
        help="Metadata column header giving cluster IDs [optional]"),
    make_option(c("-o", "--outdir"), type="character", default='.',
        help="Output directory [default %default]"),
    make_option(c("-r", "--suppress_prerotate"), action="store_true", default=FALSE,
        help="Suppress pre-rotation for optimal correlation with metadata after transformation; only relevant if metadata is supplied [default %default]"),
    make_option(c('-v','--verbose'), action='store_true', default=FALSE,
        help="Print information about execution [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), 
    args=commandArgs(trailing=TRUE))

# create output directory if needed
if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)

# LOAD DATA
pc <- read.qiime.table(opts$pcoa)
pc <- as.matrix(pc[1:(nrow(pc)-2),]) # drop last two rows of QIIME pcoa
# drop NAs in PC
keepix <- !is.na(pc[,1]) & !is.na(pc[,2])
pc <- pc[keepix,]

# load metadata
if(!is.null(opts$mapfile)){
    map <- load.qiime.mapping.file(opts$mapfile)
    qd <- remove.nonoverlapping.samples(map=map, taxa=pc)
    map <- qd$map
    pc <- qd$taxa
    # error checking
    if(!is.element(opts$category, colnames(map))){
        stop(sprintf('Gradient variable %s not found in metadata table headers\n',
                opts$category))
    }
    category <- as.numeric(map[[opts$category]])
    # drop NA values
    keep.ix <- !is.na(category)
    category <- category[keep.ix]
    pc <- pc[keep.ix,]
    map <- map[keep.ix,]
} else {
    category <- NULL
}


# DETREND
if(is.null(category) || opts$suppress_prerotate){
    resq <- get.spline.coords(pc)
} else {
    resq <- get.spline.coords.prerotate(pc, category)
}
xcoords <- resq$coords


# PLOT BEFORE/AFTER
plot.detrending(pc, xcoords, category, outdir=opts$outdir)

# REPORT CORRELATIONS BEFORE/AFTER if gradient category is present
if(!is.null(category)) write.summary.file(pc,xcoords,category,outdir=opts$outdir)

# SAVE TRANSFORMED COORDS
sink(sprintf('%s/detrended_pcoa.txt',opts$outdir))
colnames(xcoords) <- 1:ncol(xcoords)
rownames(xcoords) <- rownames(pc)
cat('pc vector number\t')
write.table(xcoords, quote=F, sep='\t')
sink(NULL)