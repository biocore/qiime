args <- commandArgs(trailingOnly=TRUE)
if(!is.element('--source_dir', args)){
    stop("\n\nPlease use '--source_dir' to specify the R source code directory.\n\n")
}
sourcedir <- args[which(args == '--source_dir') + 1]
source(sprintf('%s/util.r',sourcedir))

load.library('optparse')
load.library('DESeq2', bioconductor=TRUE)
load.library('biom')
load.library('metagenomeSeq', bioconductor=TRUE)

# make option list and parse command line
option_list <- list(
    make_option(c("--source_dir"), type="character",
        help="Path to R source directory [required]."),
    make_option(c("-i", "--input_path"), type="character",
        help="Input otu table [required]."),
    make_option(c("-o", "--out_path"), type="character", default='.',
        help="Output directory [default %default]"),
    make_option(c("-z", "--DESeq_negatives_to_zero"), type="character", default=NULL,
        help="set the negatives that result from DESeq transformation to zero")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$input_path)) stop('Please supply an otu table.')


"DESeq2"<- function(input_path, out_path, DESeq_negatives_to_zero=NULL) {
                	foo = read_biom(input_path)
                	x = as(biom_data(foo), "matrix")
                	# avoid zeros
                	x = x + 1
                	sampleTable <- data.frame(sampleName=colnames(x))
                	#Add mock design: these should not influence normalization - just required for DESeqDataSetFromMatrix
          			dds <- DESeqDataSetFromMatrix(x, sampleTable, design=~1)
                	vsmat = assay(varianceStabilizingTransformation(dds))
                	if (!is.null(DESeq_negatives_to_zero)) {
                    	vsmat[vsmat<0]=0
               		}
               		vsmat = newMRexperiment(vsmat)
               		colnames(vsmat) <- c(colnames(x))
                	write_biom(MRexperiment2biom(vsmat), out_path)
            	}

DESeq2(opts$input_path, opts$out_path, opts$DESeq_negatives_to_zero)
