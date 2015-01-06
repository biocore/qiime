args <- commandArgs(trailingOnly=TRUE)
if(!is.element('--source_dir', args)){
    stop("\n\nPlease use '--source_dir' to specify the R source code directory.\n\n")
}
sourcedir <- args[which(args == '--source_dir') + 1]
source(sprintf('%s/util.r',sourcedir))

load.library('optparse')
load.library('metagenomeSeq', bioconductor=TRUE)
load.library('biom')

# make option list and parse command line
option_list <- list(
    make_option(c("--source_dir"), type="character",
        help="Path to R source directory [required]."),
    make_option(c("-i", "--input_path"), type="character",
        help="Input otu table [required]."),
    make_option(c("-o", "--out_path"), type="character", default='.',
        help="Output directory [default %default]"),
    make_option(c("-s", "--output_CSS_statistics"), type="character", default=NULL,
        help="output CSS normalization statistics")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$input_path)) stop('Please supply an otu table.')


"CSS" <- function(input_path, out_path, output_CSS_statistics=NULL) {
			obj = load_biom(input_path)
            p = cumNormStatFast(obj)
            obj = cumNorm(obj, p = p)
            if (!is.null(output_CSS_statistics)) {
                exportStats(obj, p=p, file = file.path(output_CSS_statistics))
            }
            write_biom(MRexperiment2biom(obj, norm=TRUE, log=TRUE), out_path)
        }

CSS(opts$input_path, opts$out_path, opts$output_CSS_statistics)
