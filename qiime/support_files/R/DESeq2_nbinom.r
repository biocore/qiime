args <- commandArgs(trailingOnly=TRUE)
if(!is.element('--source_dir', args)){
    stop("\n\nPlease use '--source_dir' to specify the R source code directory.\n\n")
}
sourcedir <- args[which(args == '--source_dir') + 1]
source(sprintf('%s/util.r',sourcedir))

load.library('optparse')
load.library('DESeq2', bioconductor=TRUE)
load.library('biom')

# make option list and parse command line
option_list <- list(
    make_option(c("--source_dir"), type="character",
        help="Path to R source directory [required]."),
    make_option(c("-i", "--input_path"), type="character",
        help="Input otu table [required]."),
    make_option(c("-c", "--mapping_category"), type="character",
        help="Metadata column header giving cluster IDs [required]"),
    make_option(c("-x", "--subcategory_1"), type="character",
        help="mapping file subcategory_1, e.g. L_palm"),
    make_option(c("-y","--subcategory_2"), type="character",
        help="mapping file subcategory_2, e.g. Tongue"),
    make_option(c("-o", "--out_path"), type="character", default='.',
        help="Output directory [default %default]"),
    make_option(c("-d", "--DESeq2_diagnostic_plots"), type="character", default=NULL,
        help="provide diagnostic plots"),
    make_option(c("-e", "--outfile_diagnostic"), type="character", default=NULL,
        help="provide outfile of diagnostic plots")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$mapping_category)) stop('Please supply a mapping file header.')
if(is.null(opts$input_path)) stop('Please supply an otu table.')
if(is.null(opts$subcategory_1)) stop('Please supply a subcategory.')
if(is.null(opts$subcategory_2)) stop('Please supply a second subcategory.')


"DESeq2_nbinom" <- function(input_path, out_path, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots, outfile_diagnostic) {
            foo = read_biom(input_path)
            # remove temporary json file immediately to reduce user confusion
            file.remove(input_path, full.names=TRUE)
            x = as(biom_data(foo), "matrix")
            # avoid zeros
            x = x + 1
            colData <- data.frame(sample_metadata(foo))
            dds <- DESeqDataSetFromMatrix(x, colData, design = as.formula(paste("~",mapping_category)))
            suppressWarnings(dds <- try(DESeq(dds, quiet = TRUE), silent = TRUE))
            if (inherits(dds, "try-error")) {
                # If the parametric fit failed, try the local.
                suppressWarnings(dds <- try(DESeq(dds, fitType = "local", quiet = TRUE),
                silent = TRUE))
            if (inherits(dds, "try-error")) {
                # If local fails, try the mean
                suppressWarnings(dds <- try(DESeq(dds, fitType = "mean", quiet = TRUE),
                    silent = TRUE))
                }
            if (inherits(dds, "try-error")) {
                # If still bad, quit with error.
                return(NULL)
                }
            }
            res <- results(dds, contrast = c(mapping_category, subcategory_1, subcategory_2))
            resOrdered <- res[order(res$padj),]
            df1 <- data.frame(resOrdered)
            df1 <- cbind(OTU = rownames(df1), df1)
            write.table(df1, out_path, sep="\t", quote=F, row.names=F)
            # # #add independent filtering?
            if (!is.null(DESeq2_diagnostic_plots)) {
	            pdf(sprintf("%s", outfile_diagnostic))
	            plotMA(res, ylim = c(-3,3))
	            plotDispEsts(dds, ylim = c(1e-6, 1e1))
	            dev.off()
            }
        }

DESeq2_nbinom(opts$input_path, opts$out_path, opts$mapping_category, opts$subcategory_1, opts$subcategory_2, opts$DESeq2_diagnostic_plots, opts$outfile_diagnostic)
