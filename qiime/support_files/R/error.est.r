# File error.est.r
#
# __author__ = "Dan Knights"
# __copyright__ = "Copyright 2010, The QIIME project"
# __credits__ = ["Dan Knights"]
# __license__ = "GPL"
# __version__ = "1.2.0"
# __maintainer__ = "Dan Knights"
# __email__ = "daniel.knights@colorado.edu"
# __status__ = "Release"
 
# run with: R --vanilla --slave --args otus.txt map.txt Individual < /bio/../code/r/error.est.r

# parse arg list
argv <- commandArgs(trailingOnly=T)
source.dir <- argv[1]
x.fp <- argv[2]
map.fp <- argv[3]
categ <- argv[4]
output.dir <- argv[5]
model.names <- strsplit(argv[6],',')[[1]]

# if there are seven arguments, last one is params file: parse it
params <- NULL
if(length(argv)==7) source(argv[7])

# generate a random seed if one wasn't provided; set seed
if(!is.element('seed', names(params))) params$seed=floor(runif(1,0,1e9))
set.seed(params$seed)

# load helper files from qiime source directory
source(sprintf('%s/ml.wrappers.r',source.dir))
source(sprintf('%s/util.r',source.dir))
model.fcns <- list('random_forest'=train.rf.wrapper)

# load data
x <- t(read.table(x.fp,sep='\t',check.names=FALSE))
map <- read.table(map.fp, sep='\t',check.names=FALSE,comment.char='')
y <- as.factor(map[,categ])
names(y) <- rownames(map)

# keep only rows shared between map and data, make sure order is consistent
shared.rows <- intersect(rownames(x), names(y))
x <- x[shared.rows,]
y <- y[shared.rows]

# drop features present in < 1% of samples or in < 3 samples
rare.cutoff.by.num = 0
rare.cutoff.by.percent = ceiling(.1 * nrow(x))
if(is.element('min.num.samples',names(params))){
    rare.cutoff.by.num <- params$min.num.samples
}
if(is.element('min.percent.samples', names(params))){
    rare.cutoff.by.percent <- ceiling(params$min.percent.samples * nrow(x))
}
rare.cutoff <- max(rare.cutoff.by.num, rare.cutoff.by.percent)
x <- x[,apply(x,2,function(x) sum(x>0)) >= rare.cutoff]


# do learning, save results
for(model.name in model.names){
    subdir <- paste(output.dir,model.name,sep='/')
    load.libraries(model.name)
    model.fcn <- model.fcns[[model.name]]

    res <- tune.model(model.fcn,x,y,params)
    # save results
    save.results(res,subdir,model.name, params$seed)
}

# quit without saving history and workspace
q(runLast=FALSE)
