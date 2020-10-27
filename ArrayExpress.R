### R code from vignette source 'ArrayExpress.Rnw'

###################################################
### code chunk number 1: queryAE
###################################################
library("ArrayExpress")
sets = queryAE(keywords = "pneumonia", species = "homo+sapiens")


###################################################
### code chunk number 2: ArrayExpress-raw
###################################################
rawset = ArrayExpress("E-MEXP-1422")


###################################################
### code chunk number 3: ArrayExpress-columnsneeded
###################################################
eset = try(ArrayExpress("E-MEXP-1870"))


###################################################
### code chunk number 4: ArrayExpress-withcolumns
###################################################
eset = ArrayExpress("E-MEXP-1870",
  dataCols=list(R="ScanArray Express:F650 Mean",
    G="ScanArray Express:F550 Mean",
    Rb="ScanArray Express:B650 Mean",
    Gb="ScanArray Express:B550 Mean"))


###################################################
### code chunk number 5: getAE-full
###################################################
mexp1422 = getAE("E-MEXP-1422", type = "full")


###################################################
### code chunk number 6: ae2bioc-full
###################################################
rawset= ae2bioc(mageFiles = mexp1422)


###################################################
### code chunk number 7: getcolproc
###################################################
cn = getcolproc(mexp1422)
show(cn)


###################################################
### code chunk number 8: procset
###################################################
proset = procset(mexp1422, cn[2])


###################################################
### code chunk number 9: import (eval = FALSE)
###################################################
## AEset = ArrayExpress("E-MEXP-1416")


###################################################
### code chunk number 10: norm (eval = FALSE)
###################################################
## library("affy")
## AEsetnorm = rma(AEset)


###################################################
### code chunk number 11: fac (eval = FALSE)
###################################################
## fac = grep("Factor.Value",colnames(pData(AEsetnorm)), value=T)


###################################################
### code chunk number 12: qanorm (eval = FALSE)
###################################################
## if (suppressWarnings(require("arrayQualityMetrics", quietly=TRUE))) {
##   qanorm = arrayQualityMetrics(AEsetnorm, 
##     outdir = "QAnorm", 
##     intgroup = fac)
## }


###################################################
### code chunk number 13: limma (eval = FALSE)
###################################################
## library("limma")
## facs =  pData(AEsetnorm)[,fac]
## facs[facs[,2]=="female",2]="F"
## facs[facs[,2]=="male",2]="M"
## facs[facs[,1]=="Parkinson disease",1]="parkinson"
## facs = paste(facs[,1],facs[,2], sep=".")
## f = factor(facs)
## design = model.matrix(~0+f)
## colnames(design) = levels(f)
## fit = lmFit(AEsetnorm, design)
## cont.matrix = makeContrasts(normal.FvsM = normal.F-normal.M,             
##     parkinson.FvsM = parkinson.F-parkinson.M,  
##     Diff=(parkinson.F-parkinson.M)-(normal.F-normal.M),
##     levels=design)
## fit2 = contrasts.fit(fit, cont.matrix)
## fit2 = eBayes(fit2)
## res = topTable(fit2, coef = "parkinson.FvsM", adjust = "BH")


