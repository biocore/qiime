# John M. Gaspar (jsh58@unh.edu)
# June 2013

# This script will allow one to visualize the
#   denoising "misses" from FlowClus using
#   the lattice function "levelplot".

# output from FlowClus (-d option)
file <- "misses.csv"

library(lattice)

# load data
rep <- read.csv(file, header = FALSE)
repm <- data.matrix(rep)

# trim matrix
for (i in length(repm[1,]):1) {
  if (any(repm[,i] != 0) || any(repm[i,] != 0))
    break
}
# round up to nearest integer
while (i %% 100 != 0 && i < length(repm[1,])) {
  i <- i + 1
}
repm <- repm[1:i,1:i]

# create level boundaries -- based on max value
z <- max(repm)
reg <- c(-0.1, 0.1)
for (j in 1:7) {
  reg <- c(reg, z^(j/8))  # logarithmic: eighth powers
}
reg <- c(reg, z + 0.1)

# make arrays for plotting
int <- seq(from = 0, to = i, by = 100)
f <- log(reg[3], base = 8)
lab <- seq(from = 0, to = f * (length(reg) - 1), by = f)

# make colors using RColorBrewer
color <- RColorBrewer::brewer.pal(length(reg) - 1, "Blues")

# make levelplot
lp <- levelplot(repm,
        main = "Denoising misses",
        xlab = "Reference (cluster/node) flow value",
        ylab = "Query (read) flow value",
        scales = list(
          x = list(labels = int / 100, at = int, limits = c(0, i)),
          y = list(labels = int / 100, at = int, limits = c(0, i)),
          tck = c(1, 0)
        ),
        at = reg,
        col.regions = color,
        colorkey = list(
          col = color,
          labels = list(
            labels = c(0, 1, round(reg[3:length(reg)])),
            at = lab
          ),
          at = lab
        ),
      )
plot(lp)

# to print to a file
#jpeg("misses.jpg")
#plot(lp)
#dev.off()
