# wrapper for get.spline.coords. Pre-rotates and takes the rotation with the best
# ratio of pc1 correlation-to metadata over pc2 correlation-to-metadata
#
# nsteps is the number of rotations to make, default is 72 or 5 degrees per rotation
#
# value:
# list of
# coords: the spline coord and residual for each point
# splinef: the spline function
# proj.coords: the pc coordinates of the projection of points onto the spline
# orig.coords: the input pc coordinates
# rotated.coords: the rotated pc coordinates
"get.spline.coords.prerotate" <- function(pc, category, train.ix=NULL, nsteps=72, verbose=FALSE, ...){

    if(is.null(train.ix)){
        train.ix <- 1:nrow(pc)
    }

    thetas <- seq(0, 2 * pi, length.out = nsteps + 1)[1:(nsteps)]
    scores <- numeric(nsteps)
    best.theta <- 0
    best.resq <- NULL
    best.score <- 0
    for(i in 1:nsteps){
        theta <- thetas[i]
        
        rcoords <- rotate.2D(pc,angle=theta)
        resq <- get.spline.coords(rcoords, verbose=verbose, ...)
        xcoords <- resq$coords


        # get midpoint of pc1
        mdpt <- median(xcoords[train.ix,1])
        below <- which(xcoords[train.ix,1] <= mdpt)
        xcor.pc1 <- abs(cor(xcoords[train.ix,][,1], 
                        category[train.ix], method='spearman'))
        xcor.pc2 <- (abs(cor(xcoords[train.ix,][below,2], 
                            category[train.ix][below], method='spear')) + 
                     abs(cor(xcoords[train.ix,][-below,2], 
                            category[train.ix][-below], method='spear')))/2
        scores[i] <- xcor.pc1 / xcor.pc2
        if(scores[i] > best.score){
            best.score <- scores[i]
            best.resq <- resq
            best.theta <- theta
        }
        if(verbose){
            cat(sprintf('i = %d, theta = %.2f, score = %.3f\n', i, theta, scores[i]))
        }
    }
    return(c(best.resq,
            list(scores=scores, best.theta=best.theta, best.score=best.score)))
}


# returns the spline coordinates for one spline transformation
#
# if train.ix is not NULL, only these points are used for spline fitting
#
# value:
# list of
# coords: the spline coord and residual for each point
# splinef: the spline function
# proj.coords: the pc coordinates of the projection of points onto the spline
# orig.coords: the input pc coordinates
"get.spline.coords" <- function(pc, train.ix=NULL, ntimes=20, type=c('spline','quadratic')[2],
                        df=3, doplot=FALSE, cols='black', pch=1, verbose=FALSE){
    if(is.null(train.ix)){
        train.ix <- 1:nrow(pc)
    }

    # recursively iterate if necessary
    orig.pc <- pc
    all.m2 <- NULL
    if(ntimes > 1){
        for(i in 1:(ntimes-1)){
            res <- get.spline.coords(pc, train.ix=train.ix, type=type, df=df,
                    ntimes=1, doplot=doplot, cols=cols, pch=pch, verbose=verbose)
            all.m2 <- c(all.m2,res$m2)
            if(i==1) {
                orig.splinef <- res$splinef
                orig.proj.coords <- res$proj.coords
            }
            pc <- res$coords
        }
    }
    
    pc1 <- pc[,1]
    pc2 <- pc[,2]

    # fit detrending function to data
    x <- seq(min(pc1)-1,max(pc1)+1,.01) # fine-grained pc1 values
    n <- nrow(pc)
    # fit a spline
    if(type=='spline'){
        ss <- smooth.spline(pc1[train.ix],pc2[train.ix],df=df)
        # functions for spline
        f <- function(x) predict(ss,x)$y
        # function for spline derivative
        fx <- function(x) predict(ss,x,deriv=1)$y
    } else if(type=='quadratic'){
        mm <- lm(pc2[train.ix]~pc1[train.ix]+I(pc1[train.ix]^2))
        # functions for quadratic
        f <- function(x) {
          return(mm$coef %*% rbind(1,x,x^2))
        } 
        # function for quadratic derivative
        fx <- function(x) {
            return(mm$coef[2] + 2*mm$coef*x)
        }
    }
    
    # function for distance of point to spline at given x value
    df.pt <- function(x,func,pt) sqrt((pt[1] - x)**2 + (pt[2] - func(x))**2)
    # function for derivative of spline arc length
    fdist <- function(x) sqrt(1 + fx(x)**2)


    # for each point, find it's nearest point on the spline
    splinex <- numeric(n)
    spline.residual <- numeric(n)
    offset <- .5*diff(range(x))
    interval=c(min(x)-offset,max(x)+offset)
    for(i in 1:n){
        res <- optimize(df.pt,interval=interval,func=f,pt=pc[i,])
        splinex[i] <- res$minimum

        # if point was below the spline, make the residual negative
        if(pc[i,2] < f(splinex[i])) {
            spline.residual[i] <- -res$objective
        } else {
            spline.residual[i] <- res$objective
        }
    }

    # if this is the final iteration, store the original function and coords
    if(ntimes == 1) {
        orig.splinef <- f
        orig.proj.coords <- cbind(splinex, as.numeric(f(splinex)))
    }
    
    # traverse the points in order from left to right based on 
    # their spline x-coord, integrating distance function
    # the distance _along_ the spline is the new pc1 coord
    # the residual is the new pc2 coord
    sortixs <- sort(splinex,index=T)$ix
    spline.coord <- numeric(n) # all zeros
    for(i in 2:n){
        ix <- sortixs[i]
        previx <- sortixs[i-1]
        # try integrating along the spline
        result = tryCatch({
          dspline <- integrate(fdist, splinex[previx], splinex[ix], subdivisions=1e4)$value
        }, error = function(e) {
          dspline <- splinex[ix] - splinex[previx]
        })        
        spline.coord[ix] <- spline.coord[previx] + dspline
    }
    coords <- cbind(spline.coord,spline.residual)
    rownames(coords) <- rownames(pc)
    colnames(coords) <- c('Projected coordinates', 'Residuals')

    # calculate average change in abscissa coords
    dists <- NULL
    dists <- mean(abs((pc[,1] - mean(pc[,1])) - (coords[,1]-mean(coords[,1]))) / diff(range(pc[,1])))
    
    # m2 <- sqrt(mean((pc-coords)**2))
    m2 <- mean(dists)
    if(verbose) cat(sprintf('M^2 = %.30f\n', m2))
    res <- list(coords=coords, splinef=orig.splinef, proj.coords=orig.proj.coords, orig.coords=orig.pc, m2=c(all.m2,m2))
    if(doplot) plot.spline.coords(res, cols=cols, pch=pch)
    return(res)
}

# Writes summary of the quality of the association of the known gradient with
# coords before and after detrending
"write.summary.file" <- function(pc, xpc, category, filename='summary.txt', outdir='.'){
    pc1.cor.result <- cor.test(pc[,1], category, method='spearman', exact=FALSE)
    pc1.xcor.result <- cor.test(xpc[,1], category, method='spearman', exact=FALSE)
    
    # get midpoint of pc1
    mdpt <- median(xpc[,1])
    below <- which(xpc[,1] <= mdpt)
    pc2.cor.estimate <- (abs(cor(pc[below,2], category[below], method='spear')) + 
                          abs(cor(pc[-below,2], category[-below], method='spear')))/2
    pc2.cor.p.value.max <- max(cor.test(pc[below,2], category[below], method='spear', exact=FALSE)$p.value,
                         cor.test(pc[-below,2], category[-below], method='spear', exact=FALSE)$p.value)
    pc2.xcor.estimate <- (abs(cor(xpc[below,2], category[below], method='spear')) + 
                           abs(cor(xpc[-below,2], category[-below], method='spear')))/2
    pc2.xcor.p.value.max <- max(cor.test(xpc[below,2], category[below], method='spear', exact=FALSE)$p.value,
                              cor.test(xpc[-below,2], category[-below], method='spear', exact=FALSE)$p.value)
    # correlation of distances before vs. after
    dcor.result <- cor.test(as.numeric(dist(pc)), as.numeric(dist(xpc)),method='spear', exact=FALSE)

    sink(sprintf('%s/%s',outdir, filename))
    cat('Spearman correlation with observed gradient:\n')
    cat(sprintf('PC1 before detrending: %.4f (p=%.2e), after detrending: %.4f (p=%.2e)\n',
          -pc1.cor.result$estimate, pc1.cor.result$p.value, 
          -pc1.xcor.result$estimate, pc1.xcor.result$p.value))
    cat(sprintf('PC2 (piecewise) before detrending: %.4f (max p=%.2e), after detrending: %.4f (max p=%.3f)\n',
        pc2.cor.estimate, pc2.cor.p.value.max,
        pc2.xcor.estimate, pc2.xcor.p.value.max))    
    cat(sprintf('Spearman correlation of distances before detrending to distances after detrending: %.4f (p=%.2e).\n',
            dcor.result$estimate, dcor.result$p.value))
}

# Plots PCoA before and after detrending
"plot.detrending" <- function(pc, xpc, category=NULL, filename='PCoA_vs_projection.pdf', outdir='.'){
    filepath <- sprintf('%s/%s',outdir, filename)

    pdf(filepath,width=8,height=4)
    par(mfrow=c(1,2),mar=c(4.5,4.5,1,1), cex.lab=.8, cex.axis=.8, cex=1, lwd=1.5)
    if(is.null(category)){
        cols <- 'black'
        gradixs <- rep(1,nrow(pc))
    } else {
        cols <- brewer.pal(9,'RdYlBu') # colors for gradient
        gradixs <- get.gradient.ixs(category,9)
    }
    # plot orig coords
    plot(pc[,1], pc[,2], type='p',
            xlim=range(pc[,1]), ylim=range(pc[,2]),
            xlab='PC1',ylab='PC2', bg=cols[gradixs],pch=21)
    # plot new coords
    new.ylim <- range(pc[,2])-sum(range(pc[,2]))/2
    plot(xpc[,1], xpc[,2],type='n',col=cols[gradixs],
            xlim <- range(xpc[,1]), ylim=new.ylim, 
            xlab='Projected Coordinates',ylab='Residuals',)
    abline(0,0)
    points(xpc[,1], xpc[,2], bg=cols[gradixs], pch=21)
    dev.off()
}

# converts real values to indices between 1 and n 
"get.gradient.ixs" <- function(x,n=10){
    stdx <- x-min(x)
    stdx <- stdx/abs(diff(range(stdx)))
    gradientix <- 1 + round((n-1)*stdx)
    return(gradientix)
}

# rotates 2D coords in radians
"rotate.2D" <- function(coords, angle=2*pi/360){
    x.new = coords[,1] * cos(angle) - coords[,2] * sin(angle)
    y.new = coords[,2] * cos(angle) + coords[,1] * sin(angle)
    res <- cbind(x.new, y.new)
    rownames(res) <- rownames(coords)
    return(res)
}