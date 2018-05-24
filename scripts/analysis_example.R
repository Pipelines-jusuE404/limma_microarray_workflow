# the location of the libraries to be loaded
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("GenomicFeatures", "AnnotationDbi"))
# biocLite("ath1121501probe")
# biocLite("affyQCReport")
# biocLite("ath1121501.db")
# change the GO implementation
# biocLite("GOstats")

# library used
library(tidyverse)
library(ath1121501probe)
library(affy)
library(affyQCReport)
library(ath1121501.db)
library(ape)
library(limma)
# look for another implementation of the GO
# library(GOstats)
library(ath1121501cdf)
library(gplots)

# see the data
data_dir <- "data/CEL/"
file <- dir(data_dir)%>%
  # filter only the .cel files
  str_subset(pattern = ".cel")
################################################
# what is the identity of those cell files?
readLines(paste0(data_dir,file[1]))
# is not clear, but it contains the info relative the the array
###############################################

# ath1121501probe is a table that define the probename with sequence and position
# Print some probe sequences and their mapping positions for first Affy ID
print.data.frame(ath1121501probe[1:3,1:4])

# import all the cel file 
mydata <- ReadAffy(filenames = paste0(data_dir,file))
# what is the structure of mydata?
str(mydata)
# here are the expression data before the normalization
mydata@assayData$exprs
boxplot(mydata@assayData$exprs)
#hide the outline
boxplot(mydata@assayData$exprs,outline=F)
# data looks really not normalized
# produce a report of the data
QCReport(mydata, file="output/ExampleQC.pdf")

# Generates RMA e Background correcting Normalizing Calculating Expression
eset_rma <- rma(mydata)
# see the new distribution
boxplot(eset_rma@assayData$exprs)
# data are very well normalized

# Prints first 4 rows in data frame structure.
exprs(eset_rma)[1:4,1:2]

# exprs(eset_rma) <- 2^(exprs(eset_rma)
## Generic approach for calculating mean values for any sample combination.
# make the values as linear
mydf <- 2^exprs(eset_rma)
# what is the distribution of the result rable
boxplot(mydf,outline=F)
# remove the outline
boxplot(mydf)
# still the data are well organized

# here is an example to calculate the mean expression per probe (collapsing the replicate)
myList <- tapply(colnames(mydf), c(1,1,2,2,3,3), list)
names(myList) <- sapply(myList, paste, collapse="_")
mymean <- sapply(myList, function(x) rowMeans(mydf[,x]))
head(mymean)

## Generate MAS 5.0 P/M/A calls and combine RMA intensities, P/M/A calls
## and Wilcoxon p-values in one data frame.
eset_pma <- mas5calls(mydata)

str(eset_pma@assayData$exprs)
eset_pma@assayData$se.exprs

my_frame <- data.frame(exprs(eset_rma), exprs(eset_pma),assayDataElement(eset_pma, "se.exprs"))
# Sort columns by cel file name
my_frame <- my_frame[, sort(names(my_frame))] 
head(my_frame)

## Export results to text file that can be imported into Excel
write_csv(my_frame, "data/my_file.csv")


# Loads required annotation package.
Annot <- data.frame(ACCNUM=sapply(contents(ath1121501ACCNUM),paste,collapse=", "),
                    SYMBOL=sapply(contents(ath1121501SYMBOL),paste,collapse=", "),
                    DESC=sapply(contents(ath1121501GENENAME),paste,collapse=", "))
Annot[3:4,]
                                                                                  
## Merge annotations with expression data
all <- merge(Annot, my_frame, by.x=0, by.y=0, all=T)
all2 <- right_join(Annot%>%
             mutate(rowname=rownames(.)),
           my_frame%>%
             mutate(rowname=rownames(.))
             ,by="rowname")

## Export data to text file that can be imported into Excel
write_csv(all,"my_annot_file.csv")

# build a correlation matrix
d <- cor(2^exprs(eset_rma), method="pearson")

hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col = 4,edge.width = 2, show.node.label = T)

#####################################################################################
# analysis of the DEG

# Import sample information
targets <- readTargets("data/affy_targets.txt")
# it seems that the targets object is a dataframe
class(targets)

# load in the data notice that it is the same of the previous call at the beginning of the script
data <- ReadAffy(filenames=paste0(data_dir,targets$FileName))

# Normalization with RMA
eset <- rma(data)
                   
## If eset contains absolute intensity values like MAS5 results, then they should be transformed to log2 (or loge) values for limma. RMA/GCRMA generate log2 values and MAS5 produces absolute values.
# exprs(eset) <- log2(exprs(eset))
# Lists the analyzed file names.
pData(eset)

# Exports all affy expression values to tab delimited text file. The MAS 5.0 P/M/A calls can be retrieved with the simpleaffy package or with the affy package like this: 'eset <- mas5calls(data); write.exprs(eset, file="my_PMA.txt")'.
write.exprs(eset, file="data/affy_all.txt")

# Creates appropriate design matrix. Alternatively, such a design matrix can be created in any spreadsheet program and then imported into R.
design <- model.matrix(~ -1+factor(c(1,1,2,2,3,3)))
design
# Assigns column names.
colnames(design) <- c("group1", "group2", "group3")
# Fits a linear model for each gene based on the given series of arrays.
fit <- lmFit(eset, design)
# Creates appropriate contrast matrix to perform all pairwise comparisons. Alternatively, such a contrast matrix can be created in any spreadsheet program and then imported into R. For complex experiments one can also use this function to compute a contrast matrix with all possible pairwise comparisons.
# here are the coefficients of the contast

contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=design)

# the contrast matrix should have 1 for the row that make the contrast on the left, meaning the one that will receive the subtraction.
# the contrast matrix should have -1 for the row that make the contrast on the rigth, meaning the one that will be subtracted.
# the contrast matrix should have 0 for all the row not in the contrast

# Computes estimated coefficients and standard errors for a given set of contrasts.
fit2 <- contrasts.fit(fit, contrast.matrix)
# Computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
fit2 <- eBayes(fit2)

# Generates list of top 10 ('number=10') differentially expressed genes sorted by B-values ('sort.by=B') for each of the three comparison groups ('coef=1') in this sample set. The summary table contains the following information: logFC is the log2-fold change, the AveExpr is the average expression value accross all arrays and channels, the moderated t-statistic (t) is the logFC to its standard error, the P.Value is the associated p-value, the adj.P.Value is the p-value adjusted for multiple testing and the B-value (B) is the log-odds that a gene is differentially expressed (the-higher-the-better). Usually one wants to base gene selection on the adjusted P-value rather than the t- or B-values. More details on this can be found in the limma PDF manual (type 'limmaUsersGuide()') or on this FAQ page.

# group2-group1
topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=10)
# group3-group2
topTable(fit2, coef=2, adjust="fdr", sort.by="B", number=10)
# group3-group1
topTable(fit2, coef=3, adjust="fdr", sort.by="B", number=10)

# Exports complete limma statistics table for first comparison group ('coef=1') to tab delimited text file.
write_csv(topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000), "limma_complete.csv")

# Creates venn diagram of all changed genes with p-value equal or less than 0.05.
# it is already considering only the adj.p-value
results <- decideTests(fit2, p.value=0.05)
results
# result is a matrix of 1,0 and -1
# looking at one example ("244903_at") seems like that:
# 1 is for significant with positive FC,
# -1 is for significant with negative FC,
# 0 is for non significant.
lapply(1:3,function(x){
  topTable(fit2, coef = x, adjust = "fdr", sort.by = "B", number = 50000)%>%
    mutate(rowname=rownames(.))%>%
    filter(rowname=="244903_at")
})
vennDiagram(results)

# Filters out candidates that have P-values < 0.05 in each group ('coef=1') and provides the number of candidates for each list. These numbers should be identical with the sum of the values in each circle of the above venn diagram.
x <- topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=50000); 
y <- x[x$adj.P.Val < 0.05,];
y
# 781+787+2725+4207
# the result is the sum of all the significant genes
print("Number of genes in this list:"); nrow(y)

# Same as above but with complex filter: P-value < 0.01 AND at least 2-fold change AND expression value A > 10.
x <- topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=50000);
y <- x[x$adj.P.Val < 0.01 & (x$logFC > 1 | x$logFC < -1) & x$AveExpr > 10,];
y;
print("Number of genes in this list:"); nrow(y)

###############################################################################################################
# GO is not working this way

# This function plots heat diagram gene expression profiles for genes which are significantly differentially expressed in the primary condition (this is not a cluster analysis heat map). Genes are sorted by differential expression under the primary condition. The argument 'primary=1' selects the first contrast column in the 'results' matrix as primary condition. The plotted genes can be extracted like this 'results[results[,1]==1,]'. More information on this function can be found in the limma manual.
# results <- decideTests(fit2, p.value=0.000005);
# heatDiagram(results, fit2$coef, primary=1)
# 
# 
# affySample <- 
#   topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=100)%>%
#   mutate(ID=rownames(.))%>%
#   dplyr::select(ID)%>%
#   .[[1]]
# 
# geneSample <- na.omit(as.vector(unlist(mget(affySample, ath1121501ACCNUM, ifnotfound=NA))))
# affyUniverse <- ls(ath1121501cdf)
# geneUniverse <- as.vector(unlist(mget(affyUniverse, ath1121501ACCNUM, ifnotfound=NA)))
# 
# params <- new("GOHyperGParams", geneIds=geneSample, universeGeneIds=geneUniverse, 
#               annotation="ath1121501", ontology="MF", pvalueCutoff=0.5, conditional=FALSE, 
#               testDirection = "over")
# 
# hgOver <- hyperGTest(params)
# summary(hgOver)[c(3,7),c(1,2,5:7)]
# htmlReport(hgOver, file = "MyhyperGresult.html") 
###############################################################################################################

###############################################################################################################
# make the venn diagram with other set of filters
#The following example stores the gene identifiers for all three DEG comparisons in a list that meet the following threshold criteria: at least 2-fold change and an adjusted p-value of less then 0.01.
###############################################################################################################

# Obtain first 20 DEGs from each DEG comparison and perform hierarchical clustering on the expression matrix containing only those genes.
deglist <- sapply(1:nrow(contrast.matrix), function(x) {
  tmp <- topTable(fit2, coef=x, adjust="fdr", sort.by="B", number=Inf)%>%
    mutate(ID=rownames(.))
  affyids <- tmp[(tmp$logFC >= 1 | tmp$logFC <= -1) &
                   tmp$adj.P.Val <= 0.01, "ID"][1:20]
  })

y <- 2^exprs(eset)[unique(as.character(deglist)),]

heatmap.2(y, col=redgreen(75), scale="row", trace="none", density.info="none")

# Same as before but using correlation-based distance measures
hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete")

heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=redgreen(75),scale="row", density.info="none", trace="none")

# 1 Generate expression data with RMA, GCRMA and MAS 5.0. Create box plots for the raw data and the RMA normalized data.
# 2 Perform the DEG analysis with the limma package and determine the differentially expressed genes for each normalization data set using as cutoff an adjusted p-value of ≤0.05. Record the number of DEGs for each of the three normalization methods in a summary table.
# 3 Create for the DEG sets of the three sample comparisons a venn diagram (adjusted p-value cutoff ≤0.05).
# 4 Generate a list of genes (probe sets) that appear in all three filtered DEG sets (from B.).

sink(file = "session_info.txt");sessionInfo();sink()
