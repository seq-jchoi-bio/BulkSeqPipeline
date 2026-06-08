#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(NOISeq))

option_list <- list(
    make_option(c("-c", "--compare"), action="store_true",
                dest="compare", help="Compare conditions"),
    make_option(c("--input"), type="character", action="store",
                help="REQUIRED. File describing the input samples.",default=NULL,
                metavar="input_desc"),
    make_option(c("--info"), type="character", action="store",
                help="Optional. Table summarizing gene annotations. 
                Table must include the following: gene name, biotype, length, 
                gc, chromosome, start pos, end pos",default=NULL, metavar="file_info"),
    make_option(c("-e", "--header"), action="store_true", default = FALSE,
                dest="header", help="Do inputs have headers"),
    #make_option(c("--groups"), type="character", action="store",
    #            help="Optional. File with groups for the --info file." , default=NULL, metavar="file_groups"),
    make_option(c("-k", "--threshold"), type="integer", action="store",
                help="Optional. Threshold for the number of counts.",default=0,
                metavar="counts_threshold"),
    make_option(c("-o", "--dirOut"), type="character", action="store",
                help="Optional. Output folder.",default="./counts_qc",
                metavar="folder_output"),
    make_option("--homedir", type="character", action="store",
                help="DEVELOPMENTAL. NOT TO BE USED IN THIS VERSION", default="./",
                metavar="home_src folder")
)

opt <- parse_args(OptionParser(option_list=option_list))
HOMESRC <- opt$homedir
source(file.path(HOMESRC, "qualimapRfunctions.r"))

input.desc <- opt$input
if(is.null(input.desc)){
    stop("--input is a REQUIRED argument")
}

if(!file.exists(opt$dirOut)){
    dir.create(opt$dirOut, recursive=TRUE)
}

# cutoff for the number of counts to consider a biological feature as detected
k <- opt$threshold 

if (is.null(opt$compare)) {
    compare.conditions <- FALSE
} else {
    compare.conditions <- opt$compare
}

cm <- 1.5   # cex.main
cl <- 1.2   # cex.lab
ca <- 1     # cex.axis
cc <- 1.2   # cex

image.width <- 3*480
image.height <- 3*480
point.size <- 3*12

init.png <- function(path, width = image.width, height = image.height) {
    png(path, width, height, 
        pointsize = point.size,
        type="cairo")
}

# took it from http://gettinggeneticsdone.blogspot.de/2011/07/scatterplot-matrices-in-r.html

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

# LOAD DATA

cat("Reading input data using input description from", input.desc, "\n")
counts <- load.counts.data(input.desc, opt$header)

# Checking if data was loaded correctly
if (!all(lapply(counts, is.numeric))) {
    count_classes <- lapply(counts, class)
    print(t(count_classes))
    print(head(counts))
    stop("Input columns were not numeric. Try enabling header detection.")
}

num_samples <- ncol(counts)
cat("Loaded counts for ", nrow(counts), " features\n")
cat("Num samples:", num_samples, "\n")

expr.factors <- data.frame(Conditions = attr(counts, "factors"))

#expr.factors <- data.frame(Conditions = gl(2, num_samples/2, labels = c("C1","C2") ))
#factors <- data.frame(Tissue = c("Kidney", "Liver", "Kidney", "Liver","Kidney", "Liver","Kidney", "Liver","Kidney", "Liver"))
cat("Conditions:\n")
expr.factors
cat("\n")

cat("Compare conditions ", compare.conditions, "\n")
if ( compare.conditions == TRUE) {
    if (num_samples < 2) {
        stop("Need at least 2 samples to compare conditions")
    }
    
    if (nlevels(expr.factors[,1]) != 2)  {
        stop("Comparison can be performed only for 2 different conditions.")
    }
}  

# LOAD ANNOTATIONS

info.available <- FALSE
gene.biotypes <- NULL
gene.length <- NULL
gene.gc <- NULL
gene.loc <- NULL

if (!is.null(opt$info)){
    cat("Loading annotations from ", opt$info,"\n")
    #ann.data <- read.csv(opt$info, row.names = 1,  sep = "\t")
    ann.data <- read.csv(opt$info, sep = "\t")
    cat("Loaded annoations for ",nrow(ann.data), "genes\n")
    gene.biotypes <- ann.data[1]
    
    #cat("Gene lentgh head:", str(head(gene.length)), "\n")
    #cat(str(head(gene.biotypes)))
    gene.length <- ann.data[2]
    gene.gc <- ann.data[3]
   
    #gene.loc <- ann.data[c(4,5,6)]
    info.available <- TRUE
    
} else  {
    print("Annotation data is not available.")
}

if (info.available) {
    gene.names <- rownames(counts) 
    ann.names <- rownames(ann.data)
    intersection <- intersect(ann.names, gene.names)
    if (length(intersection) == 0) {
        warning("The gene names from counts do not correspond to the names from annotation file.")
        str(gene.names)
        str(ann.names)
        cat("Annotation based analysis is deactivated\n")
        info.available <- FALSE
        gene.biotypes <- NULL
        gene.length <- NULL
        gene.gc <- NULL
    } else  {
        #str(gene.names)
        #str(ann.names)
        cat(length(intersection),"out of",length(gene.names),"annotations from counts file found in annotations file\n")
    }
}
str(counts)
dim(counts)

cat("Init NOISeq data...\n")
ns.data <- readData(data = counts, length = gene.length, gc = gene.gc, biotype = gene.biotypes,
                    chromosome = gene.loc, factors = expr.factors )



###############################################################################

#### GLOBAL


cat("\nDraw global plots...\n")

# Counts Density

cat("Compute counts density...\n")

init.png( paste(opt$dirOut, "00_Counts_Density.png",sep="/") )
    
par(mar=c(5, 4, 4, 8) + 0.1, xpd=TRUE)

transformed.counts <- log2(counts + 0.5)

d.max <- 0
dens_data <- vector("list", num_samples)

for (i in 1:num_samples) {
    idx <- transformed.counts[,i] >= 1
    d <- density((transformed.counts[idx,i]))
    y.max <- max(d$y)
    if (y.max > d.max) {
        d.max <- y.max
    }
    dens_data[[i]] <- d
}


for (i in 1:num_samples) {
    d <- dens_data[[i]]
    if (i == 1) {
        plot( d , type='l', col = i, lwd=2, ylim=c(0,d.max),
              main="Counts density", xlab="Log2(Counts)", ylab="Density")    
    } else {
        lines(d, type='l', col = i, lwd=2)
    }
}

legend("right", title="Samples", fill=1:num_samples, legend=colnames(counts))
garbage <- dev.off()

# Correlation matrix

if (num_samples > 1) {

    cat("Compute scatterplots..\n")

    init.png( paste(opt$dirOut, "01_Scatterplot_Matrix.png",sep="/") )


    sample.size <- 5000
    if ( sample.size > nrow(transformed.counts) ) {
        sample.size <- nrow(transformed.counts)
    }

    s <- sample(1:nrow(transformed.counts), sample.size,  replace=FALSE)
    counts.sample <- transformed.counts[s,]
    pairs(counts.sample, lower.panel = panel.smooth, upper.panel = panel.cor, main = "Scatterplot Matrix")

    garbage <- dev.off()
} else {
    cat("Only one sample is provide, skipping sample comparison counts scatterplot matrix...\n")
}


# Global saturation

cat("Compute saturation..\n")
mysaturation <- dat(ns.data, k = k, ndepth = 8, type = "saturation")

init.png(paste(opt$dirOut, "02_Saturation.png", sep = "/") )
explo.plot(mysaturation, toplot = 1, samples = 1:num_samples)
garbage <- dev.off()

cat("Compute counts per biotype..\n")
counts.bio <- dat(ns.data, factor = NULL, type = "countsbio")

if (num_samples > 1) { # Workaround for 1-sample bug

# Global feature distribution plot
init.png(paste(opt$dirOut, "03_Counts_Distribution.png",sep="/"))
explo.plot(counts.bio, toplot=1, samples = NULL, plottype = "boxplot")
garbage <- dev.off()

# Global features with low count
init.png(paste(opt$dirOut, "04_Features_With_Low_Counts.png", sep="/"))
explo.plot(counts.bio, toplot=1, samples = NULL, plottype = "barplot")
garbage <- dev.off()

}

# TODO: should we include also global estimators for selected groups?

###############################################################################
#### PER SAMPLE ANALYSIS

cat("Draw per sample plots...\n\n")

if (info.available) {
    cat("Compute bio detection per sample ...\n")
    bio.detection <- dat(ns.data, k = k, type = "biodetection", factor = NULL)
    cat("Compute length bias per sample ...\n")
    length.bias <- dat(ns.data, factor = NULL, type = "lengthbias")
    cat("Compute GC bias per sample ...\n")
    gc.bias <- dat(ns.data, factor = NULL, type = "GCbias")
}

for (i in 1:num_samples) {
    
    sample.name <- colnames(counts)[i]
    cat("Processing sample",sample.name,"\n")
    
    sample.outDir <- paste(opt$dirOut, sample.name, sep = "/")
    if(!file.exists(sample.outDir)){
        dir.create(sample.outDir, recursive=TRUE)
    }

    # Saturation
    
    init.png(paste(sample.outDir,"01_Saturation.png",sep="/"))
    explo.plot(mysaturation, toplot = 1, samples = i)
    garbage <- dev.off()
    
    if (info.available) {
        
        init.png(paste(sample.outDir, "02_Bio_Detection.png", sep="/"))
        explo.plot(bio.detection, samples = i)
        garbage <- dev.off()
        
        if (num_samples > 1) { # workaround form 1-sample bug
            init.png(paste(sample.outDir, "03_Counts_Per_Biotype.png",sep="/"))
            explo.plot(counts.bio, toplot=1, samples = i, plottype = "boxplot")
            garbage <- dev.off()
        }
        
        init.png(paste(sample.outDir, "04_Length_Bias.png", sep="/"))
        explo.plot(length.bias, samples = i, toplot = 1)
        garbage <- dev.off()
        
        init.png(paste(sample.outDir, "05_GC_Bias.png", sep="/"))
        explo.plot(gc.bias, samples = i, toplot = 1)
        garbage <- dev.off()
        
    }
    
}


###############################################################################

#### PER CONDITION ANALYSIS


if (compare.conditions == TRUE) {
    
    cat ("\nPerforming comparison for conditions\n")
    
    cmp.outDir <- paste(opt$dirOut, "Comparison", sep = "/")
    if(!file.exists(cmp.outDir)){
        dir.create(cmp.outDir, recursive=TRUE)
    }
    
    cat("Compute counts per biotype for conditions..\n")
    counts.bio <- dat(ns.data, factor = "Conditions", type = "countsbio")

    init.png(paste(cmp.outDir, "01_Counts_Distribution.png",sep="/"))
    explo.plot(counts.bio, toplot=1, samples = NULL, plottype = "boxplot")
    garbage <- dev.off()
        
    init.png(paste(cmp.outDir, "02_Features_With_Low_Counts.png", sep="/"))
    explo.plot(counts.bio, toplot=1, samples = NULL, plottype = "barplot")
    garbage <- dev.off()
    
    if (info.available) {
        
        cat("Compute bio detection for conditions ...\n")
        bio.detection <- dat(ns.data, k = k, type = "biodetection", factor = "Conditions")
        
        init.png(paste(cmp.outDir, "03_Bio_Detection.png", sep="/"), width = 2*image.width)
        par(mfrow = c(1,2))
        explo.plot(bio.detection, samples = c(1,2))
        garbage <- dev.off()
        
        #init.png(paste(sample.outDir, "03_Counts_Per_Biotype.png",sep="/"))
        #explo.plot(counts.bio, toplot=1, samples = i, plottype = "boxplot")
        #garbage <- dev.off()
        
        cat("Compute length bias for conditions ...\n")
        length.bias <- dat(ns.data, factor = "Conditions", type = "lengthbias")
        
        init.png(paste(cmp.outDir, "04_Length_Bias.png", sep="/"), width = 2*image.width)
        explo.plot(length.bias, samples = NULL, toplot = 1)
        garbage <- dev.off()
        
        cat("Compute GC bias for conditions ...\n")
        gc.bias <- dat(ns.data, factor = "Conditions", type = "GCbias")    
        init.png(paste(cmp.outDir, "05_GC_Bias.png", sep="/"), width = 2*image.width)
        explo.plot(gc.bias, samples = NULL, toplot = 1)
        garbage <- dev.off()
        
    }
    
    
# #### CORRELATION PLOT
# 
# 
# if (!is.null(misdatos2)) {
#     
#     png(paste(opt$dirOut, "correlation_plot.png", sep = ""),
#         width = 3*480, height = 3*480, pointsize = 3*12)
#     
#     cor.plot.2D(misdatos1, misdatos2, noplot = 0.001, log.scale = TRUE, 
#                 xlab = paste("log2(", nom1, "+1)", sep = ""),
#                 ylab = paste("log2(", nom2, "+1)", sep = ""))
#     
#     garbage <- dev.off()
#     
# }
# 

}
    
cat("CountsQC is successfully finished!\n")









