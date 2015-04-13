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
    make_option(c("-c", "--mapping_category"), type="character",
        help="Metadata column header giving cluster IDs [required]"),
    make_option(c("-x", "--subcategory_1"), type="character",
        help="mapping file subcategory_1, e.g. Palm"),
    make_option(c("-y","--subcategory_2"), type="character",
        help="mapping file subcategory_2, e.g. Tongue"),
    make_option(c("-o", "--out_path"), type="character", default='.',
        help="Output directory [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$mapping_category)) stop('Please supply a mapping file header.')
if(is.null(opts$input_path)) stop('Please supply an otu table.')
if(is.null(opts$subcategory_1)) stop('Please supply a subcategory.')
if(is.null(opts$subcategory_2)) stop('Please supply a second subcategory.')

"fitZIG" <- function(input_path, out_path, mapping_category, subcategory_1, subcategory_2) {
            foo = read_biom(input_path)
            MGS = biom2MRexperiment(foo)
            MGS = cumNorm(MGS,p = cumNormStat(MGS))
            samplesToKeep = which(pData(MGS)[,mapping_category]%in%c(subcategory_1,subcategory_2))
            MGS = MGS[,samplesToKeep]
            MGS_category = pData(MGS)[,mapping_category]
            mod = model.matrix(~MGS_category)
            settings = zigControl(maxit=1, verbose=FALSE)
            fit = fitZig(obj=MGS, mod=mod, control=settings)
            if (is.null(observation_metadata(foo))) {
            	MRfulltable(fit, number = nrow(assayData(MGS)$counts), group=3, file=out_path)
            } else {
				res = MRfulltable(fit, number = nrow(assayData(MGS)$counts), group=3)
				sigotus = rownames(res)
				fullTaxonomyData = fData(MGS)[sigotus,]
				fullTaxonomyData = sapply(1:nrow(fullTaxonomyData), function(i){paste(format(fullTaxonomyData[i,]),collapse="; ")})
				res2 = cbind(OTU = sigotus, res,taxonomy = (fullTaxonomyData))
				write.table(res2, out_path, sep="\t", quote=F, row.names=F) 
			}
        }
        
fitZIG(opts$input_path, opts$out_path, opts$mapping_category, opts$subcategory_1, opts$subcategory_2)
