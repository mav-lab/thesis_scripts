# Take data from command line:
#####
# Infile defaults to 'Rpileup.txt'
# Outfile defaults to 'hist.txt'
# binSize defaults to 1
# hmin defaults to the lowest value in infile
# hmax defaults to the highest value in infile
#####

# Functions
getnumber <- function(prompt, default){
  cat(prompt,' (', default, ')\n')
  result <- scan(file="stdin", what="numeric", nlines=1, quiet=TRUE)
  if (length(result)>0){result<-as.numeric(result)}
  else{result <- default}
  if(!is.numeric(result) || is.na(result)) {result <- default}
  result
}

# Read command line 
args <- c(commandArgs(TRUE))
infile <- args[1]
if (is.na(infile)) {infile <- 'Rpileup.txt'}
outfile <- args[2]
if (is.na(outfile)) {outfile <- 'hist.txt'}
binSize <- as.numeric(args[3])
loLim   <- as.numeric(args[4])
hiLim   <- as.numeric(args[5])

# Read file
x <- scan(infile, quiet=TRUE);
ztemp <- c(x)
hh <- max(ztemp)
hl <- min(ztemp)

# Show info
info <- paste('Mean: ', mean(ztemp), '\n', 'Median: ', 
               median(ztemp), '\n', 'Min: ', hl, '\n', 
              'Max: ', hh, '\n')
cat('Mean: ', mean(ztemp), '\n')
cat('Median: ', median(ztemp), '\n')
cat('Min: ', hl, '\n')
cat('Max: ', hh, '\n')

# Get missing info
if (is.na(binSize)){
  binSize <- getnumber('Bin size?', 1)
}
if (is.na(loLim)){
  loLim <- getnumber('Histogram lower limit?', hl)
}
if (is.na(hiLim)){
  hiLim <- getnumber('Histogram higher limit?', hh)
}
#cat(binSize,', ',loLim,', ',hiLim,'\n')
# Calculate hist features

hlf <- binSize / 2
hmin <- loLim-hlf
hmax <- hiLim+hlf
bins <- seq(hmin,hmax,by=binSize)
if (max(bins) < hmax){bins<-c(bins[1:length(bins)-1],hmax)}
col1 <- seq(loLim,hiLim,by=binSize)
z <- ztemp[ztemp >= hmin & ztemp <= hmax]

# Save histogram
outpdft <- paste(outfile, '.jpg')
outpdf <- gsub("( )", "", outpdft)
#pdf(file=outpdf)
jpeg(filename = outpdf)
h <- hist(z, main="", breaks=bins, xlim=c(loLim,hiLim))
#plot(h, main="")
#title(main = "Histogram", sub = info, xlab = NULL, ylab = NULL, line = NA, outer = FALSE)
mytable <- matrix(c(col1,h$counts),ncol=2,byrow=FALSE)
write.table(mytable, file=outfile, row.names=FALSE, sep="\t")
