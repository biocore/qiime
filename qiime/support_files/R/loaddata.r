# reads a QIIME otu/metadata/taxon/distance table.
# Support legacy formats, where
# the header may or may not start with '#', 
# and comment lines can be anywhere in the file.
# return value is a matrix unless as.data.frame is TRUE
"read.qiime.table" <- function(filepath, as.data.frame=FALSE){
    header.index <- get.header.index(filepath)
    # read the header
    f <- file(filepath,'r')
    header <- scan(filepath, what='character', sep='\t',comment='',skip=header.index-1,quote='"',
                    nlines=1,quiet=TRUE)
    close(f)
    # read the rest of the table
    datatable <- read.table(filepath,sep='\t',skip=header.index, comment='#',quote='"',
                        head=F,row.names=1,check=FALSE,strip.white=TRUE)
    
    # set column names using header
    colnames(datatable) <- header[-1]
    
    if(!as.data.frame) datatable <- as.matrix(datatable)
    return(datatable)
}


"load.qiime.mapping.file" <- function(filepath){
    return(read.qiime.table(filepath, as.data.frame=TRUE))
}

"load.qiime.otu.table" <- function(filepath,include.lineages=FALSE){
    otus <- read.qiime.table(filepath, as.data.frame=TRUE)

    # drop "Consensus Lineage" column if present
    if(otu.table.has.metadata(colnames(otus))){
        C <- ncol(otus)
        lineages <- as.character(otus[,C])
        otus <- otus[,-C]
    } else {
        lineages <- NULL
    }
    otus <- as.matrix(t(otus))
    
    if(include.lineages){
        return(list(otus=otus,lineages=lineages))
    } else {
        return(otus=otus)
    }
}

# TRUE if last column is "Consensus Lineage" or "OTU Metadata"
"otu.table.has.metadata" <- function(headers){
    C <- length(headers)
    has.metadata <- grepl('consensus[ ]lineage|otu[ ]*metadata',
                          headers[C], ignore.case=TRUE)
    return(has.metadata)
}

# returns the index of the header line
# note: lines after the header may be comments with '#'
# read.table should therefore be called with (skip=header.index, comment='#')
"get.header.index" <- function(filepath){
    ncolumns.per.line <- NULL
    
    # read lines until the first line without a '#'
    # for each line, obtain the number of tab-delimited columns
    linecount <- 0
    start.character <- '#'
    while(start.character == '#'){
        linecount <- linecount + 1
        f <- file(filepath,'r') # open file in read mode
        line <- scan(f,what='character',skip=linecount-1,nlines=1, sep='\t',
                     quote='"', quiet=TRUE)
        close(f)
        # ncolumns is the number of entries in this line
        # not including trailing empties
        ncolumns <- max(which(sapply(line,nchar) > 0))
        ncolumns.per.line <- c(ncolumns.per.line, ncolumns)
        start.character <- substring(line[1],1,1)
    }
    
    # first non-comment line gives the number of columns
    C <- ncolumns.per.line[linecount]
    if(linecount == 1){
        # if there are no comment lines, then the first line is the header
        header.index <- 1
    } else {
        if(any(ncolumns.per.line[-linecount] == C)){
            # if there is a comment line with the correct number of columns,
            # it is the header
            header.index <- max(which(ncolumns.per.line[-linecount] == C))
        } else {
            # if there is no comment line with the correct number of columns,
            # the first non-comment line is the header
            header.index <- linecount
        }
    }

    return(header.index)
}

"load.qiime.taxon.table" <- function(filepath){
    taxa <- as.matrix(t(read.table(filepath,sep='\t',head=T,row.names=1,check=FALSE,quote='"')))
    return(taxa)
}

"load.qiime.distance.matrix" <- function(filepath){
    d <- as.matrix(read.table(filepath,sep='\t',head=T,row.names=1,check=FALSE,quote='"'))
    return(d)
}

# ensure map, data table, etc., contain the same samples in the same order
"remove.nonoverlapping.samples" <- function(map=NULL,otus=NULL,taxa=NULL,distmat=NULL){
    IDs <- NULL
    objects <- list(map=map,otus=otus,taxa=taxa,distmat=distmat)

    # find overlapping samples in all tables
    for(obj in objects){
        if(!is.null(obj)) {
            if(is.null(IDs)){
                IDs <- rownames(obj)
            } else {
                IDs <- intersect(rownames(obj), IDs)
            }
        }
    }
    
    # drop non-overlapping samples 
    for(i in 1:length(objects)){
        if(!is.null(objects[[i]])) {
            objects[[i]] <- objects[[i]][IDs,,drop=F]
            # for mapping file, drop any empty levels from factors that might
            # have occurred due to dropped samples
            if(i == 1) objects[[i]] <- droplevels(objects[[i]])
            # for distance matrix, get subset of columns too
            if(i == 4) objects[[i]] <- objects[[i]][,IDs]
        }
    }
    
    return(objects)
}
