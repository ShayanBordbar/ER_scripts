# Analysis of Bead-summary Data using beadarray
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("beadarrayExampleData", "illuminaHumanv3.db"))

library("beadarray")
library("illuminaHumanv4.db")
library("hexbin")
library(limma)
require(beadarrayExampleData)
data(exampleSummaryData)
exampleSummaryData

# getGEO gets the GEO object
GSE37386.exp.soft <- getGEO("GSE37386", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE37386.exp.soft <- GSE37386.exp.soft$GSE37386_series_matrix.txt.gz


# sampleNames gets the name of samples
sampleNames(GSE37386.exp.soft)

# featureNames gets the names of probes (rows)
featureNames(GSE37386.exp.soft)

# exprs function gets the expession matrix : probe-specific average intensities
exprs(exampleSummaryData)[1:5,1:5]
exprs(GSE37386.exp.soft)[1:8,1:5]

###

# se.exprs function gets the variability of expession : probe-specific variability
se.exprs(GSE37386.exp.soft)[1:5,1:5]

# fData function provides more information on feautures (rows)
head(fData(exampleSummaryData))
head(fData(GSE37386.exp.soft))

# pData function provides more data on samples (columns)
pData(exampleSummaryData)
pData(GSE37386.exp.soft)

#
channelNames(exampleSummaryData)
exampleSummaryData.log2 <- channel(exampleSummaryData, "G")

channelNames(GSE37386.exp.soft) # not defined here
exprs(GSE37386.exp.soft) <- log2(exprs(GSE37386.exp.soft))
###

#can subset the data using the charachter format of the names of features
randIDs <-  sample(featureNames(exampleSummaryData), 1000)
randIDs2 <- sample(featureNames(GSE37386.exp.soft), 1000)

exampleSummaryData[randIDs,]

GSE37386.exp.soft[randIDs2,]

#boxplot options: note that this is adopted from ggplot2 package, hence modifications shoukd follow a ggplot2 object
boxplot(exampleSummaryData.log2[randIDs,])
boxplot(exampleSummaryData.log2[randIDs,], what="nObservations")
boxplot(exampleSummaryData.log2[randIDs,], SampleGroup="SampleFac")
boxplot(exampleSummaryData.log2[randIDs,], probeFactor = "Status")



boxplot(GSE37386.exp.soft[1:10,])

### annotation
annotation(exampleSummaryData)
exampleSummaryData.log2 <- addFeatureData(exampleSummaryData.log2,
                                          toAdd = c("SYMBOL", "PROBEQUALITY", "CODINGZONE", "PROBESEQUENCE", "GENOMICLOCATION"))

#boxplot all probes with the gene name "ALB" as an example
ids <- which(fData(exampleSummaryData.log2)[,"SYMBOL"] == "ALB")
boxplot(exampleSummaryData.log2[ids,],
        SampleGroup = "SampleFac", probeFactor = "IlluminaID")


annotation(GSE37386.exp.soft) # annotation of my dataset is: "GPL10558". in order to work with beadarray functions I need to convert that to ExpressionSetIllumina    
#convert to ExpressionSetIllumina
summaryData <- as(GSE37386.exp.soft, "ExpressionSetIllumina")
head(fData(summaryData))
# mark regular vs negatie control probes
fData(summaryData)$Status <-
  ifelse(fData(summaryData)$PROBEQUALITY=="No match","negative","regular" )
fData(summaryData)$Status[1:20]

#I don't know what exactly this is but it says:a detection p-value (Detection Pval) estimates the probability of a gene being detected above the background level
Detection(summaryData) <- calculateDetection(summaryData, status=fData(summaryData)$Status)

#This is the normalization function which can be done using different methods, my data is already log2 transformed and quantile normalized
summaryData.norm <- normaliseIllumina(summaryData,method="quantile",
                                      status=fData(summaryData)$Status, transform = "log2")

# MA plot: to see if normalization is rquired. most points should lie around y=0 line
# M=\log _{2}(R/G)= \log _{2}(R)-\log _{2}(G)
# A={\frac  12}\log _{2}(RG)={\frac  12}(\log _{2}(R)+\log _{2}(G))
mas <- plotMA(summaryData,do.log=FALSE)
mas
##Added lines on the y axis
mas + geom_hline(yintercept=c(-1.5,1.5),col="red",lty=2)
##Added a smoothed line to each plot
mas+ geom_smooth(col="red")
##Changing the color scale
mas + scale_fill_gradient2(low="yellow",mid="orange",high="red")

# FIlTERATION of bad data
# We recommend removing probes assigned a ‘Bad’ or ‘No match’ quality score after normalization

ids <- as.character(featureNames(summaryData))
#getting the quality of probes from the corresponding annotation database. in this case illuminaHumanv4
qual <- unlist(mget(ids, illuminaHumanv4PROBEQUALITY, ifnotfound=NA))
table(qual)

#remove probes with bad , nomatch or NA quality
rem <- qual == "No match" | qual == "Bad" | is.na(qual)
summaryData.filt <- summaryData[!rem,]
dim(summaryData.filt)

#### differential expression
#example data
rna <- factor(pData(exampleSummaryData)[,"SampleFac"])
design <- model.matrix(~0+rna)
colnames(design) <- levels(rna)
aw <- arrayWeights(exprs(exampleSummaryData.log2), design)
aw
fit <- lmFit(exprs(exampleSummaryData.log2), design, weights=aw)
contrasts <- makeContrasts(UHRR-Brain, levels=design)
contr.fit <- eBayes(contrasts.fit(fit, contrasts))
topTable(contr.fit, coef=1)

#my data
# want to look at differential expression between (knock out + E2) and (knock out + vehicle) in 6 or 24 hours

aa <- paste(pData(summaryData.filt)[,42],pData(summaryData.filt)[,43], pData(summaryData.filt)[,44] ,sep = "_")
summaryData.filt.ordered <- summaryData.filt[, sort(aa, index.return = T)$ix]
aa2 <- paste(pData(summaryData.filt.ordered)[,42],pData(summaryData.filt.ordered)[,43], pData(summaryData.filt.ordered)[,44] ,sep = "_")

aafac <- factor(aa2)
aadesign <-  model.matrix(~0+aafac)

colnames(aadesign) <- levels(aafac)
aw <- arrayWeights(exprs(summaryData.filt.ordered), aadesign)
aw
fit <- lmFit(exprs(summaryData.filt.ordered), aadesign, weights=aw)
contrasts <- makeContrasts(siGreb_24hrs_E2 - siGreb_control_vehicle, levels=aadesign)
contr.fit <- eBayes(contrasts.fit(fit, contrasts))
topTable(contr.fit, coef=1, number = 100)

contrasts <- makeContrasts(siNT_6hrs_E2 - siNT_control_vehicle, levels=aadesign)
contr.fit <- eBayes(contrasts.fit(fit, contrasts))
topTable(contr.fit, coef=1, number = 100)

rownames(topTable(contr.fit, coef=1, number = 100))
ada <- fData(summaryData.filt.ordered)[rownames(topTable(contr.fit, coef=1, number = 100)),"Symbol"]

length(unique(ada))




