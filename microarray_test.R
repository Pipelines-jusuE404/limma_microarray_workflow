#######################
## Homework R Script ##
#######################
## How to run this script:
# 	source("homework_script.R")

## Dependencies:
# 	R working directory needs to contain
# 	all *.cel files and targets.txt file
#       R packages:
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("affy", "limma", "gcrma")
         
         ## RMA Commands
         ## (A) Normalization: RMA
         library(affy); library(limma); library(gcrma) # Loads required libraries.
         targets <- readTargets("targets.txt") # Import targets information.
         data <- ReadAffy(filenames=targets$FileName) # Import expression raw data and stores them as AffyBatch object.
         probedata <- data.frame(ID=probeNames(data), PM=pm(data), MM=mm(data)) # Access proble level data.         
         eset_rma <- rma(data) # Normalizes the data with 'rma' function and assigns them to exprSet object.
         # exprs(eset) <- log2(exprs(eset)) # Only MAS5 stores absolute intensities. GCRMA and RMA methods store log2 intensities.
         # pData(eset_rma) # Lists the analyzed file names.
         # write.exprs(eset_rma, file="rma_all.xls") # Exports all affy expression values to tab delimited text file. 
         ## Create Box Plots for Raw Data and Normalized Data
         pdf(file="raw_boxplot.pdf"); boxplot(data, col="red", main="Raw Data"); dev.off() # Generates boxplot for un-normalized log intensity values.
         pdf(file="rma_boxplot.pdf"); boxplot(data.frame(exprs(eset_rma)), col="blue", main="RMA Normalized Data"); dev.off() # Generates boxplot for RMA normalized log intensity values.
         
         ## (B) DEG Analysis for RMA Data
         design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2,3,3,3))) # Creates appropriate design matrix.
         colnames(design) <- c("S1", "S2", "S3") # Assigns nicer column names.
         contrast.matrix <- makeContrasts(S2-S1, S3-S2, S3-S1, levels=design) # Creates appropriate contrast matrix for pairwise comparisons.
         fit <- lmFit(eset_rma, design) # Fits a linear model for each gene based on the given series of arrays.
         fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
         fit2 <- eBayes(fit2) # Computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
         rma_deg_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
         rma_deg_result <- rma_deg_result[rma_deg_result$adj.P.Val<=0.05,]
         write.table(rma_deg_result, "rma_deg_result.xls", col.names = NA, quote=FALSE, sep="\t")
         
         ## (C) Create Venn Diagram for RMA Data
         rma_venn <- decideTests(fit2, p.value=0.05)
         pdf(file="rma_venn.pdf"); vennDiagram(rma_venn); dev.off() # Creates venn diagram of all changed genes with p-value equal or less than 0.05.
         
         ## (A) Normalization: GCRMA
         eset_gcrma <- gcrma(data) # Normalizes the data with 'gcrma' function and assigns them to exprSet object.
         
         ## (B) DEG Analysis for GCRMA Data
         fit <- lmFit(eset_gcrma, design) # Fits a linear model for each gene based on the given series of arrays.
         fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
         fit2 <- eBayes(fit2) # Computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
         gcrma_deg_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
         gcrma_deg_result <- gcrma_deg_result[gcrma_deg_result$adj.P.Val<=0.05,]
         write.table(gcrma_deg_result, "gcrma_deg_result.xls", quote=FALSE, col.names = NA, sep="\t")
         
         ## (C) Create Venn Diagram for GCRMA Data
         gcrma_venn <- decideTests(fit2, p.value=0.05)
         pdf(file="gcrma_venn.pdf"); vennDiagram(gcrma_venn); dev.off() # Creates venn diagram of all changed genes with p-value equal or less than 0.05.
         
         ## (A) Normalization: MAS5
         eset_mas5 <- mas5(data) # Normalizes the data with 'mas5' function and assigns them to exprSet object.
         exprs(eset_mas5) <- log2(exprs(eset_mas5)) # Only MAS5 stores absolute intensities. Other methods store log2 intensities.
         
         ## (B) DEG Analysis for MAS5 Data
         fit <- lmFit(eset_mas5, design) # Fits a linear model for each gene based on the given series of arrays.
         fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
         fit2 <- eBayes(fit2) # Computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
         mas5_deg_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
         mas5_deg_result <- mas5_deg_result[mas5_deg_result$adj.P.Val<=0.05,]
         write.table(mas5_deg_result, "mas5_deg_result.xls", quote=FALSE, col.names = NA, sep="\t")
         
         ## (C) Create Venn Diagram for MAS5 Data
         mas5_venn <- decideTests(fit2, p.value=0.05)
         pdf(file="mas5_venn.pdf"); vennDiagram(mas5_venn); dev.off() # Creates venn diagram of all changed genes with p-value equal or less than 0.05.
         
         ## (D) Identifiy the Overlap Between the Three Methods
         overlap2 <-(merge(rma_deg_result, gcrma_deg_result, by.x = "ID", by.y = "ID", all = FALSE))
         overlap3 <-(merge(overlap2, mas5_deg_result, by.x = "ID", by.y = "ID", all = FALSE))
         write.table(overlap3, "overlap.xls", quote=FALSE, col.names = NA, sep="\t")
         