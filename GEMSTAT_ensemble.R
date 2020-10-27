################LIBRARIES
library(gplots)
library(RColorBrewer)
########################################################################################################################
########################################################################################################################
########################################                              ##################################################
########################################          FUNCTIONS           ##################################################
########################################                              ##################################################
########################################################################################################################
########################################################################################################################
#1. EnsembleParWriter 2. GEMSTAT_input_constructor_Ensemble 3. GEMSTAT_output_Reader_ensemble 4. getMode 5. plot_ensemble_output
#6. LogFileReader 7. PerformanceFilter 8. EnhWeightReader 9. RegElemWeightHeatMap 10. RegElemWeightCoordinate
#11. ParameterHeatMap 12. PerformanceHeattMap 13. CompareModelPlot 14. logSimilarityCheck 15. FilterForModels 16. PredictedvsRealCor
# 17. PlotPercentGenes 18. PurturbationInputCreator 19. Purturbation_output_Reader

########################################################################################################################
########################################################################################################################

EnsembleParWriter <- function(TFnames, par.RangeMat, nu.sample, .coopertingTFIndex){
  # Create an ensemble of initial Par files for GEMSTAT
  # inputs to this function:
  # 1. TFnames : a charachter vector containing the name of TFs
  # 2. par.RangeMat : a matrix containing the range of each of the parameter and if the parameter is log scale or not: 
  # each row indicates parameter, first column is lower bound, second column is upper bound, third column is 1 if the ranges are in log scale and 0 otherwise
  # 3. nu.sample : number of samples per compartment: we will have 2^number of parameters compartments, and "nu.sample" indicates how many samples shuld be retrieved from a compartment
  # 4. .coopertingTFIndex is a matrix where nrow = number of interacting pairs, ncol= 2 . Each row has the index of interacting TFs in one interaction (index in TF expression matrix)
  
  number2binary = function(number, noBits) {
    # function creates the binary version of a number
    binary_vector = rev(as.numeric(intToBits(number)))
    if(missing(noBits)) {
      return(binary_vector)
    } else {
      binary_vector[-(1:(length(binary_vector) - noBits))]
    }
  }
  
  range<- par.RangeMat
  Par.nu= nrow(range);
  M<- matrix(0, nrow=2^Par.nu, ncol=Par.nu);
  for(i in 0:(2^Par.nu-1))
  {
    M[i+1,] = number2binary(i,Par.nu)
  }
  
  compartmentRanges <- list()
  for(i in 1:nrow(M)){
    rangeFile <- matrix(0, nrow=Par.nu, ncol=3);
    rangeFile[,3] <- range[, 3];
    for(j in 1:Par.nu){
      if(M[i, j]==0){
        rangeFile[j, 1] <- range[j, 1]
        rangeFile[j, 2] <- (range[j, 2]+range[j, 1])/2
      }
      else if(M[i,j]==1){
        rangeFile[j, 1] <- (range[j, 2]+range[j, 1])/2
        rangeFile[j, 2] <- range[j, 2]
      }
      else{
        rangeFile[j, 1] <- range[j, 1]
        rangeFile[j, 2] <- range[j, 2]
      }
    }
    compartmentRanges[[i]] <- rangeFile
  }
  
  #now sample from the ranges for each compartment
  compartmentSamples <- list()
  for(comp in 1:length(compartmentRanges)){
    ranges <- compartmentRanges[[comp]]
    Nsample <- nu.sample
    n <- nrow(ranges)
    sample <- matrix(0, nrow=Nsample, ncol=n)
    for(i in 1:n){
      if(ranges[i, 3]==0){
        sample[, i] <- runif(Nsample, ranges[i, 1], ranges[i, 2])
      }else{
        sample[, i] <- 10^(runif(Nsample, ranges[i, 1], ranges[i, 2]))
      }
    }
    compartmentSamples[[comp]] <- sample
  }
  #write parameter par for each
  for(comp in 1:length(compartmentRanges)){
    for(samp in 1:nu.sample){
      curSample <- compartmentSamples[[comp]][samp, ]
      for(TFnu in 1:length(TFnames)){
        #write TF parameters
        cat(paste(TFnames[TFnu], curSample[(TFnu)*2 -1],
                  curSample[(TFnu)*2], sep = "\t"),
            file = paste0("par", comp, "_", samp,".txt"),
            sep = "\n", append = T)
      }
      #write basal
      cat(paste0("basal_transcription = ", curSample[(TFnu)*2 + 1]), file = paste0("par", comp, "_",samp,".txt"), sep = "\n", append = T)
      #write coop
      for(cptf in 1:nrow(.coopertingTFIndex)){
        cat(paste(TFnames[.coopertingTFIndex[cptf,1]], TFnames[.coopertingTFIndex[cptf,2]], curSample[length(TFnames)*2 + 1 + cptf], sep = "\t"), file = paste0("par", comp, "_", samp,".txt"), sep = "\n", append = T)
      }
    }
  }
}
#####################################
#####################################
#example
aa <- matrix(0L, nrow = 6, ncol = 3)
aa[1,] <- c(-2, 4, 1)
aa[2,] <- c(-2, 4, 1)
aa[3,] <- c(-2, 4, 1)
aa[4,] <- c(-2, 4, 1)
aa[5,] <- c(-2, 3, 1)
aa[6,] <- c(-2, 3, 1)
EnsembleParWriter(TFnames = c("ESR1", "RARA"), par.RangeMat = aa, nu.sample = 20, .coopertingTFIndex = rbind(c(1,1),c(1,2)))
###############################################################################################################
###############################################################################################################

###############################################################################################################
###############################################################################################################
GEMSTAT_input_constructor_Ensemble <- function(clusterIndex,
                                               gene_expression_Mat, 
                                               .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                               .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                               .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                               .promoter.sequence=Promoter.Sequence.Char,
                                               singleEnhancer=F,
                                               Manual_enhancer_Seq=character(0),
                                               .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                               .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                                               .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                               experiment.nu,
                                               .par.exp.nu = experiment.nu,
                                               .motifs=TF.motifs.count,
                                               .TF_expression=TF.Expression.Quantile.For.GEMSTAT,
                                               .max.enh.length = numeric(0),
                                               ..a=character(0), ..o="DIRECT", coopertingTFIndex=numeric(0),
                                               ..i=character(0), ..r=character(0), ..oo=character(0), ..mc=character(0),
                                               ..p=character(0), ..rt=character(0), ..na=character(0), ..ct=character(0), ..sigma=character(0)
                                               ,.par.RangeMat, .nu.sample, writeParameters = T){
  # get the cluster and create all gemstat inputs and copies them to veda, also creates the job file
  
  #This function gets as inout:
  # 1) clusterIndex <- a list where each entry is index of genes in one cluster
  # 2) gene_expression_Mat <- matrix of gene expression which the indexes in clusterIndex refer to. It has to be named
  #  entries for : Function to write the bag of regulatory elements of each gene
  # 3).enhancer.per.gene : a list containing list of enhancers for each gene eg: Vicinity100kb.Enhancer.plusPromoter.By.gene , for all genes: It will be subsetted using gene names
  # 4) .Enhancer.Sequence : matirx of all enhancers: MCFEnhancers
  # 5).promoter.promoter.int : a list containing list of promoters interacting with each gene eg: Promoter.gene.PromoterIntList: for all genes: It will be subsetted using gene names
  # 6) .promoter.sequence : sequence of the promoters of genes in the same order as input1 and input 3--->this should be named by genes, 
  # IT IS IMPORTANT THAT THE Third AND Fifth AND Forth ENTRIES ARE NAMED, (with corresponding genes names)
  # IT IS IMPORTANT THAT THE Third AND Fifth  ENTRIES SHOULD HAVE THE SAME LENGTH, WHICH IS THE NUMBER OF GENES OF INTEREST
  # singleEnhancer is True if we want to work with one enhancer per gene and we specify the index of that enhancer for each gene in Manual_enhancer_Index
  # Manual_enhancer_Seq is an character vector containing the sequence of enhancers in input: .Enhancer.Sequence. One enhancer per gene, ordered by genes in expression matrix. only read when singleEnhancer is True
  # 7) dir: the directory in which the files will be stored
  # 8) experiment.nu is the number of experiment for which this data is being created
  # 8_2 ) .par.exp.nu : is the number of experiment from which the parameters are being used, default set to current experiment
  # 9) .motifs : motifs in the format to be written 
  # 10).TF_expression : TF experession matrix 
  # 11) .max.enh.length : maximum length of considered enhancers
  # 12) last row are gemstat job parameters used in Hal_job_writer:
  # ..a annFile
  # ..o modelOption
  # coopertingTFIndex is a matrix where nrow = number of interacting pairs, ncol= 2 . Each row has the index of interacting TFs in one interaction (index in TF expression matrix)
  # ..i factorInfoFile
  # ..r repressionFile 
  # ..oo objOption 
  # ..mc maxContact 
  # ..p parFile # if you want to use parameter file it has to be created and its address be added here
  # ..rt repressionDistThr 
  # ..na nAlternations 
  # ..ct coopDistThr 
  # ..sigma factorIntSigma
  
  # .par.RangeMat parameters for creating ensemble parameter files
  # .nu.sample parameters for creating ensemble parameter files
  # writeParameters : boolean, if True creates and writes parameters, if False it doesn't
  
  # This function creates: for each cluster: Input and Output folder, populates the Input folder with sequence, gene and TF expression and motif data.
  # This function creates a bash file for copying the created directory to veda
  # This function creates a job file to be run in hal
  # This function uses five other functions: Hal_job_writer, WriteFastaOfBag, ExpressionWriter, MotifWriter, bash_directory
  
  prev_dir <- getwd()
  setwd(dir = .dir)
  dir.create(paste0(.dir,"/Experiment_", experiment.nu))
  for(clus in 1:length(clusterIndex)){
    print(paste0("Creating cluster ", clus, " input files ..."))
    # creating directory for this cluster
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus))
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus, "/Inputs"))
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus, "/Outputs"))
    # writing the expression matrix of this cluster
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus, "/Inputs", "/Gene_expression"))
    setwd(paste0(.dir,"/Experiment_", experiment.nu,
                 "/Cluster_", clus, "/Inputs", "/Gene_expression"))
    print("Gene Expression")
    current_exp_mat <- gene_expression_Mat[clusterIndex[[clus]],]
    ExpressionWriter(current_exp_mat,
                     output.File.Name = paste("ExpressionMat_",
                                              as.character(clus), sep=""))
    # write the bag of regulatory elements
    print("Sequence")
    setwd(.dir)
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus, "/Inputs", "/Sequence"))
    setwd(paste0(.dir,"/Experiment_", experiment.nu,
                 "/Cluster_", clus, "/Inputs", "/Sequence"))
    Index_in_orig <- match(rownames(current_exp_mat), names(.enhancer.per.gene))
    if(singleEnhancer){#if inputting single enhancers per gene
      for(reg.el in 1:length(Manual_enhancer_Seq)){
        #Manual_enhancer_Index is an integer vector containing the index of enhancers in input: .Enhancer.Sequence. One enhancer per gene, ordered by genes in expression matrix
        write.table(c(paste(">", paste(rownames(current_exp_mat)[reg.el], "1", sep=":"), sep=""),
                      Manual_enhancer_Seq[reg.el]),
                    file=paste0(paste("RegulatoryElementsSequence", as.character(clus), sep="_"), ".fa"),
                    sep="\n", row.names=F, col.names=F, quote=F, append=T)
      }
      
    }else{ #if inputting genes and writing sequences based on (enhacner, promoter, gene) associations
      WriteFastaOfBag(enhancer.per.gene=.enhancer.per.gene[Index_in_orig],
                      Enhancer.Sequence=.Enhancer.Sequence,
                      promoter.promoter.int=.promoter.promoter.int[Index_in_orig],
                      promoter.sequence=.promoter.sequence,
                      output.File.Name=paste("RegulatoryElementsSequence", as.character(clus), sep="_"),
                      max.enh.length=.max.enh.length)
    }

    setwd(.dir)
    
    # write the Motifs
    print("Motifs")
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus, "/Inputs", "/Motifs"))
    setwd(paste0(.dir,"/Experiment_", experiment.nu,
                 "/Cluster_", clus, "/Inputs", "/Motifs"))
    MotifWriter(motif.List=.motifs, output.File.Name = "motifs")
    
    print("TF Expression")
    # write the TF Expression
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus, "/Inputs", "/TF_Expression"))
    setwd(paste0(.dir,"/Experiment_", experiment.nu,
                 "/Cluster_", clus, "/Inputs", "/TF_Expression"))
    ExpressionWriter(.TF_expression ,output.File.Name = "TF_Expression_quantile")
    
    # write parameters
    if(writeParameters){
      print("Writing parameter files ...")
      dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                        "/Cluster_", clus, "/Inputs", "/Parameters"))
      setwd(paste0(.dir,"/Experiment_", experiment.nu,
                   "/Cluster_", clus, "/Inputs", "/Parameters"))
      EnsembleParWriter(TFnames = rownames(.TF_expression), par.RangeMat = .par.RangeMat, nu.sample = .nu.sample, .coopertingTFIndex = coopertingTFIndex)
    }else{
      print(paste("skip writing parameters, Using parameters from experiment", .par.exp.nu))
    }
  }
  # write the script for hal job
  setwd(paste0(.dir,"/Experiment_", experiment.nu))
  
  if (length(coopertingTFIndex) > 0){
    print("Writing the coop file")
    #create the cooperation file and feed to Hal_job_writer 
    TF.RN <- rownames(.TF_expression)
    for(cptf in 1:nrow(coopertingTFIndex)){
      cat(paste(TF.RN[coopertingTFIndex[cptf,1]], TF.RN[coopertingTFIndex[cptf,2]], sep = "\t"), file = "coop.txt", sep = "\n", append = T)
    }
    ..c <- paste0(.sharedmountsdir,"Experiment_", experiment.nu,"/coop.txt")
  }else{
    ..c <- character(0)
  }
  
  print("creating the hal job file ...")
  Hal_job_writer(exp.nu=experiment.nu, cluster.nu=length(clusterIndex),
                 seqName="RegulatoryElementsSequence",
                 expressionName="ExpressionMat",
                 Shared_dir= .sharedmountsdir,
                 GEMSTAT_dir = .GEMSTAT_dir,
                 par.exp.nu = .par.exp.nu,
                 Ensemble = T,
                 home_dir = paste0(.dir,"/"),
                 .a=..a , .o=..o , .c=..c, .i=..i, .r=..r, .oo=..oo, .mc=..mc, .p=..p, .rt=..rt, .na=..na, .ct=..ct, .sigma=..sigma)
  # create a .submit from the written hal job
  sys_com <- paste0(paste0(paste0(paste0("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl"," ","Experiment_",
                                         experiment.nu),"job "),
                           length(clusterIndex) * (2^nrow(.par.RangeMat)) * .nu.sample / 10),
                    paste0(paste0(paste0(paste0(" tmpjob_", experiment.nu),
                                         " >job_exp_"), experiment.nu), ".submit"))
  cat(c("#!/bin/bash", "\n"),file=paste0("hal_sub_creator_", "exp_", experiment.nu), sep="")
  cat(sys_com, file=paste0("hal_sub_creator_", "exp_", experiment.nu), sep="", append = T)
  
  sys_com <- paste0("chmod +x ", paste0(" hal_sub_creator_", "exp_", experiment.nu))
  system(sys_com)
  #write a bash file to copy the created files to shared-mounts
  bash_directory(exp.nu=experiment.nu,
                 root_dir="tabebor2@veda.cs.illinois.edu:/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble",
                 source_dir=.dir)
  #make the created bash file executable
  sys_com <- paste0(paste0("chmod +x ", paste0(paste0(.dir,"/Experiment_", experiment.nu),
                                               "/bash_directory_exp_" )), experiment.nu)
  system(sys_com)
  print("Copying to veda ...")
  sys_com <- paste0(paste0(paste0(.dir,"/Experiment_", experiment.nu),
                           "/bash_directory_exp_" ), experiment.nu)
  system(sys_com)
  # reset the working directory to where it was
  setwd(prev_dir)
}
#########################################################################################################
#########################################################################################################
#example
aa <- matrix(0L, nrow = 6, ncol = 3)
aa[1,] <- c(-2, 4, 1)
aa[2,] <- c(-2, 4, 1)
aa[3,] <- c(-2, 4, 1)
aa[4,] <- c(-2, 4, 1)
aa[5,] <- c(-2, 3, 1)
aa[6,] <- c(-2, 3, 1)
GEMSTAT_input_constructor_Ensemble(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed))),
                                               gene_expression_Mat = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed, 
                                               .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                               .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                               .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                               .promoter.sequence=Promoter.Sequence.Char,
                                               .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                               .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                                               .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                               experiment.nu = 11,
                                               .motifs=TF.motifs.TimeSeries.count,
                                               .TF_expression=TF.Expression.TimeSeries,
                                               .max.enh.length = numeric(0),
                                               ..o="DIRECT", coopertingTFIndex=rbind(c(1,1),c(1,2))
                                               , ..oo="SSE", ..ct=13, ..sigma=character(0)
                                               ,.par.RangeMat = aa, .nu.sample= 20)

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
GEMSTAT_output_Reader_ensemble <- function(.exp.nu,
                                           .number.Of.Clusters,
                                           .root_dir="tabebor2@veda.cs.illinois.edu:/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble",
                                           .dest_dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble"
                                           ,copyOutput = T){
  # this function mounts to veda (root_dir), copies the output file of the given experiments for all of 
  # its clusters to the directory of that experiment on my laptop (dest_dir)
  # It then reads the output expression file, and outputs the expressio, SSE and ### plots for each gene if needed ### plotting not impelemented
  # Inputs:
  # 1) .exp.nu : number of the experiment
  # 2) .number.Of.Clusters = number of clusters in that experiment
  # 3) .root_dir : the root directory in veda where Experiment folders are located
  # 4) .dest_dir : the destination directory on my laptop where experiment folders are located.
  ####### not impelemented 5) if you want the output plotted
  # write a bash file to copy the created files to shared-mounts
  prev_dir <- getwd()
  setwd(paste0(.dest_dir, "/Experiment_", .exp.nu))
  on.exit(setwd(prev_dir))
  if(copyOutput){
    bash_directory(exp.nu=.exp.nu,
                   root_dir=.root_dir,
                   source_dir=.dest_dir,
                   Read.output = T,
                   number.Of.Clusters=.number.Of.Clusters,
                   file_name = "bash_directory_output_")
    # make the created bash file executable
    sys_com <- paste0("chmod +x ","bash_directory_output_exp_", .exp.nu)
    system(sys_com)
    print("Copying from veda ...")
    sys_com <- paste0("./bash_directory_output_exp_", .exp.nu)
    system(sys_com)
  }
  AllResults <- list()
  tmpFile <- list.files(paste0(.dest_dir, "/Experiment_",
                               .exp.nu,"/Cluster_",1,"/Outputs/"),pattern = "*.outEns")
  par_number <- length(tmpFile)
  SSEResults <- matrix(nrow =par_number , ncol = .number.Of.Clusters)
  rownames(SSEResults) <- tmpFile
  for(clus in 1:.number.Of.Clusters){
    AllResults[[clus]] <- list()
    print(paste0("reading output for cluster_",clus,"..."))
    output_file_names <- list.files(paste0(.dest_dir, "/Experiment_",
                                           .exp.nu,"/Cluster_",clus,"/Outputs/"),pattern = "*.outEns")
    OutPut_clus <- list() #list contaning the outputs of this cluster
    for(output_file in 1:length(output_file_names)){
      OutPut_clus[[output_file]] <- read.table(paste0(.dest_dir, "/Experiment_",
                                                      .exp.nu,"/Cluster_",clus,"/Outputs/", output_file_names[output_file]),
                                               header = T, sep = "\t")
      cur_Prediction <- OutPut_clus[[output_file]][seq(2, nrow(OutPut_clus[[output_file]]), 2), 2:ncol(OutPut_clus[[output_file]])]
      cur_real <- OutPut_clus[[output_file]][seq(1, nrow(OutPut_clus[[output_file]]), 2), 2:ncol(OutPut_clus[[output_file]])]
      cur_Prediction <- do.call(cbind, cur_Prediction)
      cur_real <- do.call(cbind, cur_real)
      rownames(cur_Prediction) <- unique(OutPut_clus[[output_file]][, 1])
      rownames(cur_real) <- rownames(cur_Prediction)
      cur_SSE <- sqrt(rowSums((cur_real - cur_Prediction)^2)/ncol(cur_real))
      names(cur_SSE) <- rownames(cur_real) 
      AllResults[[clus]][[output_file]] <- list(prediction=cur_Prediction,
                                                TrueValue=cur_real,
                                                SSE = cur_SSE)
      SSEResults[output_file, clus] <- mean(cur_SSE)
    }
    names(AllResults[[clus]]) <- output_file_names
  }#end of loop over clusters
  
  
  return(list(Allresults = AllResults, averageSSE = SSEResults))
}
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
#example
Exp11.GEMSTAT.output <- GEMSTAT_output_Reader_ensemble(.exp.nu = 11, .number.Of.Clusters = 1)
par(mfrow = c(1,1), mar = c(4,4,4,4))
hist(Exp11.GEMSTAT.output$averageSSE[,1], breaks = 100)
sort(Exp11.GEMSTAT.output$averageSSE[,1],decreasing = F ,index.return = T)$ix
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
#plotting output
plot_ensemble_output <- function(ensembleOutput,
                                 ClusterNumber=1, paramnumber=0,
                                 OneGeneManyModel = F, GeneNumber = 0, Rank = 1, varAmong = integer(0), plotAmong = integer(0),
                                 groupvsgroup = F, group1 = integer(0), group2 = integer(0)){
  # ensembleOutput : is the output of GEMSTAT_output_Reader_ensemble
  # ClusterNumber : number of cluster you want to plot
  # paramnumber : number of parameter (index in second output of GEMSTAT_output_Reader_ensemble) you want to plot
  # OneGeneManyModel : if True plots the prediction of one gene in all parameter settings
  # GeneNumber : index of the gene which you want to be plotted in OneGeneManyModel setting
  # Rank : if GeneNumber is not defined, it chooses the gene that has highest variance of performance in differnet models, Rank is the rank of variance among all genes.
  # varAmong : is a integer vector containing the index of the models among which you want to calculate the variance of genes prediction
  # plotAmong : is a integer vector containing the index of the models which you want to plot
  # groupvsgroup if True tries to find the genes that are differently predicted between two groups of parameters
  # group1 : integer vector of all indecis of models in group1 
  # group2 : integer vector of all indecis of models in group2
  prev.Par <- par()
  ngenes <- nrow(ensembleOutput$Allresults[[ClusterNumber]][[1]]$prediction)
  nconditions <- ncol(ensembleOutput$Allresults[[ClusterNumber]][[1]]$prediction)
  nmodels <- length(ensembleOutput$Allresults[[ClusterNumber]])
  if(OneGeneManyModel){
    if(GeneNumber == 0){
      if(groupvsgroup){#if you want to find the genes that are most different between two groups
        # get the predicted expression of all genes in each group
        group1AllPrediction <- array(rep(0, length(group1)*ngenes*nconditions), dim = c(length(group1), ngenes, nconditions))
        group2AllPrediction <- array(rep(0, length(group2)*ngenes*nconditions), dim = c(length(group2), ngenes, nconditions))
        for(i in 1:length(group1)){
          group1AllPrediction[i,,] <- ensembleOutput$Allresults[[ClusterNumber]][[group1[i]]]$prediction
        }
        for(i in 1:length(group2)){
          group2AllPrediction[i,,] <- ensembleOutput$Allresults[[ClusterNumber]][[group2[i]]]$prediction
        }
        # get the median expression value for each condition, for each gene
        group1MedianPrediction <- apply(X = group1AllPrediction, MARGIN = c(2,3), median)
        group2MedianPrediction <- apply(X = group2AllPrediction, MARGIN = c(2,3), median)
        # compute the SSE between median of the two groups
        group12SSE <- rowSums((group1MedianPrediction - group2MedianPrediction)^2)
        # sort based on difference of medians
        group12SSERank <-sort(group12SSE, decreasing = T, index.return = T)$ix
        # choose the Rank th gene of the sorted list
        ChosenGene <- group12SSERank[Rank]
      }else{#if gene is not provided and not in groupvsgroup mode, sorts the genes based on variance and chooses the Rank th gene
        modelIndices <- c(1:length(ensembleOutput$Allresults[[ClusterNumber]]))
        if(length(varAmong) > 0){
          modelIndices <- varAmong
        }
        GeneSSE <- matrix(nrow = length(ensembleOutput$Allresults[[ClusterNumber]][[1]]$SSE), ncol = length(modelIndices))
        for (i in 1:length(modelIndices)){
          GeneSSE[,i] <- ensembleOutput$Allresults[[1]][[modelIndices[i]]]$SSE
        }
        GeneSSEVar <- apply(GeneSSE,1,var)
        GeneSSEVarRank <-sort(GeneSSEVar, decreasing = T, index.return = T)$ix
        ChosenGene <- GeneSSEVarRank[Rank]
        }
    }else{
      ChosenGene <- GeneNumber
    }
    #plot performance of the gene in differenct models
    ChosenGeneReal <- ensembleOutput$Allresults[[ClusterNumber]][[1]]$TrueValue[ChosenGene, ]
    print(paste("Chosen gene index is",ChosenGene))
    if(groupvsgroup){# plotting  two groups of models against each other for one gene
      par(mfrow = c(8, 8), mar = c(0.1, 0.1, 0.1, 0.1))
      for(curmodel in 1:length(group1)){
        #plot for the chosen gene in each selected model
        plot(ensembleOutput$Allresults[[ClusterNumber]][[group1[curmodel]]]$prediction[ChosenGene, ],
             ylim = c(0,1), col = 2,
             xaxt = "n", yaxt = "n", main ="", ylab = "", xlab = "",
             type = "l", lwd = 2)
        lines(ChosenGeneReal, ylim = c(0,1), col = 3)
        aa <- ensembleOutput$Allresults[[ClusterNumber]][[group1[curmodel]]]$SSE[ChosenGene]
        text((length(ChosenGeneReal) - 2), 0.9 , format(round(aa, 2), nsmall = 2), cex = 0.5) # write the size of the SSE
        text(4, 0.1 , group1[curmodel], cex = 0.5) # write the number of model
      }
      for(curmodel in 1:length(group2)){
        #plot for the chosen gene in each selected model
        plot(ensembleOutput$Allresults[[ClusterNumber]][[group2[curmodel]]]$prediction[ChosenGene, ],
             ylim = c(0,1), col = 4,
             xaxt = "n", yaxt = "n", main ="", ylab = "", xlab = "",
             type = "l", lwd = 2)
        lines(ChosenGeneReal, ylim = c(0,1), col = 3, bg="grey")
        aa <- ensembleOutput$Allresults[[ClusterNumber]][[group2[curmodel]]]$SSE[ChosenGene]
        text((length(ChosenGeneReal) - 2), 0.9 , format(round(aa, 2), nsmall = 2), cex = 0.5) # write the size of the SSE
        text(4, 0.1 , group2[curmodel], cex = 0.5) # write the number of model
      }
    }else{ #plotting certain models
      par(mfrow = c(11, 12), mar = c(0.1, 0.1, 0.1, 0.1))
      modelIndecisPlot <- c(1:length(ensembleOutput$Allresults[[ClusterNumber]]))
      if(length(plotAmong) > 0){
        modelIndecisPlot <- plotAmong
      }
      for(curmodel in 1:length(modelIndecisPlot)){
        #plot for the chosen gene in each selected model
        plot(ensembleOutput$Allresults[[ClusterNumber]][[modelIndecisPlot[curmodel]]]$prediction[ChosenGene, ],
             ylim = c(0,1), col = 2,
             xaxt = "n", yaxt = "n", main ="", ylab = "", xlab = "",
             type = "l", lwd = 2)
        
        lines(ChosenGeneReal, ylim = c(0,1), col = 3)
        aa <- ensembleOutput$Allresults[[ClusterNumber]][[modelIndecisPlot[curmodel]]]$SSE[ChosenGene]
        text((length(ChosenGeneReal) - 2), 0.9 , format(round(aa, 2), nsmall = 2), cex = 0.5) # write the size of the SSE
        text(4, 0.1 , modelIndecisPlot[curmodel], cex = 0.5) # write the number of model
      }
    }

  }else{ # if plotting all genes for one model: paramnumber
    # get the prediction of the model for all genes
    cur_Prediction <- ensembleOutput$Allresults[[ClusterNumber]][[paramnumber]]$prediction
    # get the tru value of expression for all genes
    cur_real <- ensembleOutput$Allresults[[ClusterNumber]][[paramnumber]]$TrueValue
    # get the SSE of each gene
    cur_SSE <- ensembleOutput$Allresults[[ClusterNumber]][[paramnumber]]$SSE
    par(mfrow = c(11, 12), mar = c(0.1, 0.1, 0.1, 0.1))
    for(gene in 1:nrow(cur_Prediction)){#plot for each gene
      plot(cur_Prediction[gene,],
           ylim = c(0,1), col = 2,
           xaxt = "n", yaxt = "n", main ="", ylab = "", xlab = "",
           type = "l", lwd = 2)
      
      lines(cur_real[gene,], ylim = c(0,1), col = 3)
      aa <- cur_SSE[gene]
      text((ncol(cur_Prediction) - 2), 0.9 , format(round(aa, 2), nsmall = 2), cex = 0.5) # write the size of the SSE
      text(2, 0.1 , gene, cex = 0.5) # write the index of the gene
      
    }
  }
  #set plot parameters back to what it was
  par(mfrow = prev.Par$mfrow, mar = prev.Par$mar)
}
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
LogFileReader <- function(logdir,
                          removeDups = F,
                          nu.TFs=2,
                          nu.coop.pairs=2,
                          parameterNames=c("ESR1_Binding","ESR1_txp","RARA_binding","RARA_txp","basal","ESR1_ESR1","RARA_ESR1")){
  # read log file and extract parameters and performance
  # input
  #1) logdir : the address of the directory where log files live
  #2) removeDups : if True removes duplicates and almost duplicates
  #3) nu.TFs is the number of TFs used in this model
  #4) nu.coop.pairs : number of cooperative pairs defined
  #5) parameterNames: a character vector containing the name of parameters in this order: 
  # c("TF1_Binding", "TF1_txp", "TF2_binding", "TF2_txp", ... , "basal", "coop1_coop1", "coop2_coop2", ...)
  
  # outputs a list of two entries, first one contains parameters, second one contains performance
  # first entry is a matrix nrow: number of models which successfully finished execution, ncol = number of parameters
  # second entry is a numeric vector, each entry is performance of a model, corresponding to the parameter in the first entry
  
  #write a bash file that gets a log file and only keeps the lines corresponding to parameters and performance
  # cat(c("#!/bin/bash", "\n"),file=paste0("logFilter_exp_", exp.nu), sep="")
  # cat(c("cd ",logdir, "\n"),file=paste0("logFilter_exp_", exp.nu), sep="", append=T)
  # cat(c(paste0("for filename in ",logdir, "/*.logEns; do"), "\n"), file=paste0("logFilter_exp_", exp.nu), sep="", append=T)
  # cat("$filename")
  
  logfilenames <- list.files(logdir, pattern = "*.Filtered")
  parMatrix <- matrix(0L, nrow = length(logfilenames), ncol = (nu.TFs*2 + nu.coop.pairs + 1))
  rownames(parMatrix) <- logfilenames
  colnames(parMatrix) <- parameterNames
    # c("ESR1_Binding","ESR1_txp","RARA_binding","RARA_txp","basal","ESR1_ESR1","RARA_ESR1")
  perfVector <- numeric(length(logfilenames))
  names(perfVector) <- logfilenames
  for(i in 1:length(logfilenames)){
    aa <- readLines(paste0(logdir,"/",logfilenames[i]))
    #reading binding and txp parameters for each TF
    for (cur_TF in 1:nu.TFs){
      parMatrix[i,(cur_TF*2 - 1):(cur_TF*2)] <- as.numeric(strsplit(aa[1 + cur_TF],split = "\t")[[1]][2:3])
    }
    #reading basal
    parMatrix[i,(2 * nu.TFs + 1)] <- as.numeric(strsplit(aa[2 + nu.TFs],split = " ")[[1]][3])
    #reading coop parameters
    for(cur_coop in 1:nu.coop.pairs){
      parMatrix[i,(2 * nu.TFs + 1 + cur_coop)] <- as.numeric(strsplit(aa[2 + nu.TFs + cur_coop],split = "\t")[[1]][3])
    }
    #reading performance calculated by GEMSTAT
    perfVector[i] <- as.numeric(strsplit(aa[2 + nu.TFs + nu.coop.pairs + 1],split = " ")[[1]][3])
  }
  parMatrix <- parMatrix[!is.na(parMatrix[,1]), ]
  perfVector <- perfVector[!is.na(perfVector)]
  if(removeDups){
    parMatrix <- parMatrix[!duplicated(parMatrix), ]
    perfVector <- perfVector[!duplicated(parMatrix)]
    #remove rows that are almost duplicates (sum of abolute differences less than 0.1)
    dup <- numeric(nrow(parMatrix))
    for(parrow in 1:(nrow(parMatrix) - 1)){
      for(j in (parrow+1):nrow(parMatrix)){
        if(dup[j]==0){
          aa <- parMatrix[parrow, ] - parMatrix[j, ]
          aaa <- sum(abs(aa)/ncol(parMatrix))
          #print(aaa)
          if(aaa < 1){
            dup[j] <- parrow
          }
        }
      }
    }
    print(paste((length(dup) - sum(dup %in% 0)), " out of ", length(dup), " remaining models marked as almost duplicated and removed."))
    parMatrix <- parMatrix[dup %in% 0, ]
    perfVector <- perfVector[dup %in% 0]
  }

  result = list(parameters = parMatrix, performance = perfVector)
  return(result)
}
####################################################################################################################################################
#example
aa <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_14/Cluster_1/Outputs")
aa <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_14/Cluster_1/Outputs")
####################################################################################################################################################
####################################################################################################################################################
logSimilarityCheck <- function(LogFileReader.Output, logscale = F, thresh = 1){
  #This function is the duplicate checking part of LogFileReader function
  parMatrix <- LogFileReader.Output$parameters
  perfVector <- LogFileReader.Output$performance
  parMatrix <- parMatrix[!duplicated(parMatrix), ]
  perfVector <- perfVector[!duplicated(parMatrix)]
  if(logscale){
    CheckMat <- log10(parMatrix)
  }else{
    CheckMat <- parMatrix
  }
  #remove rows that are almost duplicates (sum of abolute differences less than 0.1)
  dup <- numeric(nrow(CheckMat))
  for(parrow in 1:(nrow(CheckMat) - 1)){
    for(j in (parrow+1):nrow(CheckMat)){
      if(dup[j]==0){
        aa <- CheckMat[parrow, ] - CheckMat[j, ]
        aaa <- sum(abs(aa)/ncol(CheckMat))
        #print(aaa)
        if(aaa < thresh){
          dup[j] <- parrow
        }
      }
    }
  }
  print(paste((length(dup) - sum(dup %in% 0)), " out of ", length(dup), " remaining models marked as almost duplicated and removed."))
  parMatrix <- parMatrix[dup %in% 0, ]
  perfVector <- perfVector[dup %in% 0]
  result = list(parameters = parMatrix, performance = perfVector, removed = which(! dup %in% 0))
  return(result)
}
####################################################################################################################################################
####################################################################################################################################################
#example
aa <- logSimilarityCheck(Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5)
aaa <- logSimilarityCheck(Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5, logscale = T, thresh = 0.1)
####################################################################################################################################################
####################################################################################################################################################
FilterForModels <- function(ModelIndex = integer(0),GEMSTAT_output_Reader_ensemble.output, EnhWeightReader.output, 
                            basedOnSimilarityCheck = T, .LogFileReader.Output, .logscale = T, .thresh = 0.1 ){
  # This function filters the models in output of GEMSTAT_output_Reader_ensemble,EnhWeightReader and LogFileReader based on
  # the index of models given: ModelIndex, or basedOnSimilarityCheck: in which case runs similarity check of parameters and filters based on that
  # inputs: 
  # 1. ModelIndex : integer vector containing the index of models we want to keep
  # 2. GEMSTAT_output_Reader_ensemble.output : output of GEMSTAT_output_Reader_ensemble function
  # 3. EnhWeightReader.output : output of 3.EnhWeightReader function
  # 4. basedOnSimilarityCheck : if True runs logSimilarityCheck on .LogFileReader.Output and keeps the passed models
  # 5. .LogFileReader.Output :  output of .LogFileReader function
  # 6. .logscale : parameter of logSimilarityCheck function: if the comparison be run of log scale
  # 7. .thresh : parameter of logSimilarityCheck function: threshold of similarity
  
  if(basedOnSimilarityCheck){
    aaa <- logSimilarityCheck(LogFileReader.Output = .LogFileReader.Output, logscale = .logscale, thresh = .thresh)
    ChosenModels <- setdiff(c(1:nrow(.LogFileReader.Output$parameters)), aaa$removed)
    remove(aaa)
  }else{
    ChosenModels <- ModelIndex
  }
  # first filter enhancer contributions
  EnhWeightReader.output.Filtered <- EnhWeightReader.output[ChosenModels]
  # next filter output
  GEMSTAT_output_Reader_ensemble.output.Filtered <- GEMSTAT_output_Reader_ensemble.output
  GEMSTAT_output_Reader_ensemble.output.Filtered$averageSSE <- GEMSTAT_output_Reader_ensemble.output$averageSSE[ChosenModels]
  GEMSTAT_output_Reader_ensemble.output.Filtered$Allresults[[1]] <- GEMSTAT_output_Reader_ensemble.output$Allresults[[1]][ChosenModels]
  # next filter parameters
  .LogFileReader.Output.Filtered <- .LogFileReader.Output
  .LogFileReader.Output.Filtered$parameters <- .LogFileReader.Output$parameters[ChosenModels, ]
  .LogFileReader.Output.Filtered$performance <- .LogFileReader.Output$performance[ChosenModels]
  results = list(output = GEMSTAT_output_Reader_ensemble.output.Filtered, param = .LogFileReader.Output.Filtered, enh = EnhWeightReader.output.Filtered)
  return(results)
}

####################################################################################################################################################
####################################################################################################################################################
#example
aaa <- FilterForModels(ModelIndex = c(1:10),
                       GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.output.sigmoid_no_softmax_5,
                       EnhWeightReader.output = Exp17.GEMSTAT.enhWeight.sigmoid_no_softmax_5,
                       basedOnSimilarityCheck = F,
                       .LogFileReader.Output = Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5)
####################################################################################################################################################
####################################################################################################################################################
#write a function to filter a group of models based on paerformance on one or a group of genes
PerformanceFilter <- function(ensembleOutput ,cluster = 1, modelIndex=numeric(0),
                              geneIndex, perfThresh = 0.18, parameters = numeric(0),
                              sortModels = T, strictRemove = F){
  #filters models based on performance on a  set of genes. inputs:
  # 1) ensembleOutput : is the output of GEMSTAT_output_Reader_ensemble
  # 2) cluster : index of cluster to be considred
  # 3) modelIndex : integer vector containing the index of models to be considred for filtering
  # 4) index of gene or genes that filtering is based on them
  # 5) perfThresh : threshold on performance
  # 6) parameter: output of LogFileReader
  # 7) if sortModel is True: sorts all given models based on their average performance on the given group of genes
  # 8) if strictRemove = T removes the models which have worse than threshold performance on any of the given genes
  # 9) if sortModel mode: returns a matrix which first column is the index of the input models sorted based on perfomance, second column is average performance
  ensembleSize <- length( ensembleOutput$Allresults[[cluster]])
  TotalgeneNum <- length( ensembleOutput$Allresults[[cluster]][[1]]$SSE)
  if (length(modelIndex) == 0){
    modelIndex <- c(1:ensembleSize)
  }
  if(strictRemove){
    PassFilter <- numeric(length(modelIndex))
    for(curModel in 1:length(modelIndex)){
      passCond <- (ensembleOutput$Allresults[[cluster]][[modelIndex[curModel]]]$SSE[geneIndex] < perfThresh)
      if(all(passCond)){
        PassFilter[curModel] <- 1
      }else{
        print(paste("for model",modelIndex[curModel]," ",(length(passCond) -sum(passCond)), " out of ",length(passCond), " didn't pass the test"))
      }
    }
    return(modelIndex[PassFilter %in% 1])
  }else if(sortModels){
    PerformanceNotSorted <- numeric(length(modelIndex))
    for(curModel in 1:length(modelIndex)){
      PerformanceNotSorted[curModel] <- mean(ensembleOutput$Allresults[[cluster]][[modelIndex[curModel]]]$SSE[geneIndex])
    }
    PerformanceSorted <- sort(PerformanceNotSorted, decreasing = F, index.return = T)
    return(cbind(PerformanceSorted$ix,  PerformanceSorted[[1]]))
  }
}
####################################################################################################################################################
####################################################################################################################################################
#example
aa <- PerformanceFilter(ensembleOutput = Exp16.GEMSTAT.output ,cluster = 1,
                        modelIndex=numeric(0), geneIndex = c(16, 13),
                        perfThresh = 0.18)
mean(Exp16.GEMSTAT.output$Allresults[[1]][[2601]]$SSE[aa$ClusterIndex[[1]]])
####################################################################################################################################################
####################################################################################################################################################
#write a function to read the enhancer weights given the .EnhWeight file
EnhWeightReader <- function(Directory){
  # Directory is the directory in which .EnhWeight files live
  #output is a list of lists: first each list correspond to a model, then each list within the model correspond to a gene (contains a numeric vector (length: number of enhancers for that gene), each entry is the weight of that enhancer for the gene in this model)
  EnhWfilenames <- list.files(Directory,pattern = "*.EnhWeight")
  EnhWList <- list()
  
  for(i in 1:length(EnhWfilenames)){
    EnhWList[[i]] <- list()
    #read all lines of the file
    All_Lines <- readLines(paste0(Directory, "/", EnhWfilenames[i]))
    #split each line
    All_Lines_Splitted <- strsplit(All_Lines, split = "_")
    splitted_Len <- unlist(lapply(All_Lines_Splitted, length))
    #remove the lines which have length of splitted equal to one
    # All_Lines_Splitted <- All_Lines_Splitted[!(splitted_Len %in% 1)]
    #remove the lines which have length of splitted equal to two
    # All_Lines_Splitted <- All_Lines_Splitted[!(splitted_Len %in% 2)]
    #keep the lines which have length of splitted equal to six
    All_Lines_Splitted <- All_Lines_Splitted[(splitted_Len %in% 6)]
    # splitted_Len <- splitted_Len[!(splitted_Len %in% 1)]
    # splitted_Len <- splitted_Len[!(splitted_Len %in% 2)]
    splitted_Len <- splitted_Len[(splitted_Len %in% 6)]
    if(sum(splitted_Len %in% 6) != length(splitted_Len)){
      print(splitted_Len)
      print(paste(sum(splitted_Len %in% 6),"vs", length(splitted_Len)))
      print("something is wrong with splitting")
    }
    #convert the All_Lines_Splitted list into a dataframe:
    All_Lines_Splitted_df <- do.call(rbind, All_Lines_Splitted)
    All_Lines_Splitted_df <- All_Lines_Splitted_df[!(duplicated(All_Lines_Splitted_df)), ]
    
    #get the forth element of each entry which is the gene name
    unique_genes <- unique(All_Lines_Splitted_df[, 4])
    #
    for(gene in 1:length(unique_genes)){
      EnhWList[[i]][[gene]] <- as.numeric(All_Lines_Splitted_df[All_Lines_Splitted_df[, 4] %in% unique_genes[gene], 6])
    }
    names(EnhWList[[i]]) <- unique_genes
  }
  aa1 <- strsplit(EnhWfilenames, split = "\\.")
  a1 <- unlist(lapply(aa1, '[[', 2))
  names(EnhWList) <- a1
  return(EnhWList)
}
####################################################################################################################################################
####################################################################################################################################################
#example
Exp16.GEMSTAT.enhWeight.sigmoid_1_softmax_10 <- EnhWeightReader("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_16/Cluster_1/Outputs")
####################################################################################################################################################
####################################################################################################################################################
RegElemWeightHeatMap <- function(.EnhWeightReaderOutput, .col_vec=col_vector,
                                 .modelIndex=c(1:length(.EnhWeightReaderOutput)),
                                 .removeNonFunctionalRegs=F,
                                 .nonFuncThreshold=0.3,
                                 .exportplot=T,
                                 .dendrogram = "column",
                                 .Colv=T,.Rowv = F,
                                 .filename = "enhancerHeatmap.png",
                                 setmaxto1RestZero = F,
                                 .ColSideColors = character(0),
                                 .RowSideCol = character(0),
                                 enhancerIndex= integer(0)){
  # setmaxto1RestZero : if True for each gene sets the weight of the enhancer with max weight to 1 and rest of its enhancers to zero (used depending on the GEMSTAT obj func)
  # .ColSideColors is charachter vector length equal to the number of models, determines the color on the sidebar for each model.
  # enhancerIndex : integer vector containing the index of enhancers that you want to be plotted. note that every other input should be provided for the full enhancer list, then in the function it will subset to use the enhancerIndex enhancers
  library(RColorBrewer)
  library(gplots)
  if(length(.ColSideColors) == 0){
    .ColSideColors = rep("white", length(.modelIndex))
  }
  RegElNumber <- sum(unlist(lapply(.EnhWeightReaderOutput[[1]], length)))
  RegWeightMat <- matrix(nrow = RegElNumber, ncol = length(.modelIndex))
  for(model in 1:length(.modelIndex)){
    RegWeightMat[, model] <- unlist(.EnhWeightReaderOutput[[.modelIndex[model]]])
  }
  colnames(RegWeightMat) <- as.character(.modelIndex)
  rownames(RegWeightMat) <- c(1:nrow(RegWeightMat))
  ##
  if(setmaxto1RestZero){
    maxAssign1 <- function(inpvec){
      #gets a numeric or int vector, returns a vector of the same length, with all zeros except the maximum which is one
      a <- which.max(inpvec)
      aa <- inpvec * 0
      aa[a] <- 1
      return(aa)
    }
    genassoc <- integer(0)
    geneEnhLen <- unlist(lapply(.EnhWeightReaderOutput[[1]], length))
    for(gene in 1:length(.EnhWeightReaderOutput[[1]])){
      genassoc <- c(genassoc, rep(gene, geneEnhLen[gene]))
    }
    for(gene in 1:length(.EnhWeightReaderOutput[[1]])){
      curInd <- which(genassoc == gene)
      if(length(curInd) == 1){
        RegWeightMat[curInd, ] <- rep(1, ncol(RegWeightMat))
      }else{
        RegWeightMat[curInd, ] <- apply(X = RegWeightMat[curInd, ], MARGIN = 2, FUN = maxAssign1)
      }
    }
  }
  ##
  my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
  if(length(.RowSideCol) == 0){
    print("no .RowSideCol defined")
    #create side colors for each row: same color enhancers correspond to the same gene
    RowSideCol <- numeric(0)
    geneEnhLen <- unlist(lapply(.EnhWeightReaderOutput[[1]], length))
    for(gene in 1:length(.EnhWeightReaderOutput[[1]])){
      RowSideCol <- c(RowSideCol, rep(.col_vec[(gene %% length(.col_vec)) + 1], geneEnhLen[gene]))
    }
    RowSideCol <- as.character(RowSideCol)
  }else{
    print(".RowSideCol defined")
    RowSideCol <- .RowSideCol
  }

  if(.removeNonFunctionalRegs){
    rowmax <- apply(RegWeightMat, 1, max)
    funcRegIndex <- which(rowmax > .nonFuncThreshold)
    RegWeightMat <- RegWeightMat[funcRegIndex, ]
    RowSideCol <- RowSideCol[funcRegIndex]
  }
  
  #export or not?
  if(.exportplot){
    png(.filename,    # create PNG for the heat map        
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
  }
  #plot the heatmap using gplot heatmap.2
  if(length(enhancerIndex) == 0){
    enhancerIndex <- c(1:nrow(RegWeightMat))
  }
  
  heatmap.2(RegWeightMat[enhancerIndex, ],
            #cellnote = RegWeightMat,  # same data set for cell labels
            main = "Reg elements", # heat map title
            notecol="black",      # change font color of cell labels to black
            # density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            # margins =c(12,9),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            dendrogram=.dendrogram,     # only draw a col dendrogram
            RowSideColors = RowSideCol[enhancerIndex], # grouping row-variables into different categories
            ColSideColors = .ColSideColors,
            Rowv = .Rowv,
            Colv = .Colv,
            scale = "none"
  ) 
  if(.exportplot){
    dev.off() 
  }
}
####################################################################################################################################################
####################################################################################################################################################
# example
aaaa <- sort(Exp17.GEMSTAT.output.sigmoid_no_softmax_5$averageSSE, index.return = T, decreasing = F)$ix
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp17.GEMSTAT.enhWeight.sigmoid_no_softmax_5,.modelIndex = aaaa[1:20], .removeNonFunctionalRegs = F, .exportplot = T)

RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$enh,
                     .removeNonFunctionalRegs = T, .nonFuncThreshold = 0.1,
                     .exportplot = T,
                     .modelIndex = aad,
                     .dendrogram = "none",
                     .Colv = F, .Rowv = F,
                     setmaxto1RestZero = T)


####################################################################################################################################################
####################################################################################################################################################

RegElemWeightCoordinate <- function(RegulatoryCoordinateExtractorOutput,
                                    EnhWeightReaderOutput,
                                    modelIndex = c(1:length(EnhWeightReaderOutput)),
                                    plotheatmap = T, col_vec = col_vector, exportplot = T, filename = "enhancerHeatmap.png",
                                    removeNonFunctionalRegs = F, nonFuncThreshold = 0.2, setmaxto1RestZero = F){
  # RegulatoryCoordinateExtractor is in GRO-Seq.R script
  # this function gets the output of RegulatoryCoordinateExtractor function and output of EnhWeightReader function and combines them to
  # add a column to RegulatoryCoordinateExtractor output per model
  # modelIndex is the index of the model in EnhWeightReaderOutput
  # plotheatmap : if True plots a heatmap of regulatory elements contribution, each row is a regulatory element, each column is a model
  # col_vec : a vector of colors
  # exportplot : if True exports the plot to a png file
  # filename : name of the file to export to
  # removeNonFunctionalRegs : remove the regulatory elements which do not contribute to any of the models.
  # setmaxto1RestZero : if True for each gene sets the weight of the enhancer with max weight to 1 and rest of its enhancers to zero (used depending on the GEMSTAT obj func)
  
  # matrix that contains the weight of each reg element in each model
  RegWeightMat <- matrix(nrow = nrow(RegulatoryCoordinateExtractorOutput), ncol = length(modelIndex))
  for(model in 1:length(modelIndex)){
    RegWeightMat[, model] <- unlist(EnhWeightReaderOutput[[modelIndex[model]]])
  }
  colnames(RegWeightMat) <- as.character(modelIndex)
  rownames(RegWeightMat) <- c(1:nrow(RegWeightMat))
  ######
  if(setmaxto1RestZero){
    maxAssign1 <- function(inpvec){
      #gets a numeric or int vector, returns a vector of the same length, with all zeros except the maximum which is one
      a <- which.max(inpvec)
      aa <- inpvec * 0
      aa[a] <- 1
      return(aa)
    }
    genassoc <- integer(0)
    geneEnhLen <- unlist(lapply(EnhWeightReaderOutput[[1]], length))
    for(gene in 1:length(EnhWeightReaderOutput[[1]])){
      genassoc <- c(genassoc, rep(gene, geneEnhLen[gene]))
    }
    for(gene in 1:length(EnhWeightReaderOutput[[1]])){
      curInd <- which(genassoc == gene)
      if(length(curInd) == 1){
        RegWeightMat[curInd, ] <- rep(1, ncol(RegWeightMat))
      }else{
        RegWeightMat[curInd, ] <- apply(X = RegWeightMat[curInd, ], MARGIN = 2, FUN = maxAssign1)
      }    }
  }
  ######
  result <- cbind(RegulatoryCoordinateExtractorOutput, RegWeightMat)
  
  ###############################plotting############################
  
  if (plotheatmap){
    library(RColorBrewer)
    library(gplots)
    my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
    #create side colors for each row: same color enhancers correspond to the same gene
    RowSideCol <- numeric(0)
    geneEnhLen <- unlist(lapply(EnhWeightReaderOutput[[1]], length))
    for(gene in 1:length(EnhWeightReaderOutput[[1]])){
      RowSideCol <- c(RowSideCol, rep(col_vec[(gene %% length(col_vec)) + 1], geneEnhLen[gene]))
    }
    RowSideCol <- as.character(RowSideCol)
    # if removeNonFunctionalRegs option is on, modify RegulatoryCoordinateExtractorOutput, RegWeightMat and RowSideCol
    if(removeNonFunctionalRegs){
      rowmax <- apply(RegWeightMat, 1, max)
      funcRegIndex <- which(rowmax > .nonFuncThreshold)
      # funcRegIndex <- which(rowSums(RegWeightMat) > nonFuncThreshold)
      RegulatoryCoordinateExtractorOutput <- RegulatoryCoordinateExtractorOutput[funcRegIndex, ]
      RegWeightMat <- RegWeightMat[funcRegIndex, ]
      RowSideCol <- RowSideCol[funcRegIndex]
    }

    #export or not?
    if(exportplot){
      png(filename,    # create PNG for the heat map        
          width = 5*300,        # 5 x 300 pixels
          height = 5*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8)        # smaller font size
    }
    #plot the heatmap using gplot heatmap.2
    heatmap.2(RegWeightMat,
              #cellnote = RegWeightMat,  # same data set for cell labels
              main = "Reg elements", # heat map title
              notecol="black",      # change font color of cell labels to black
              # density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              # margins =c(12,9),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              #breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="column",     # only draw a col dendrogram
              RowSideColors = RowSideCol, # grouping row-variables into different categories
              Rowv = F,
              scale = "none"
              ) 
    if(exportplot){
      dev.off() 
    }
  }else{
    if(removeNonFunctionalRegs){
      rowmax <- apply(RegWeightMat, 1, max)
      funcRegIndex <- which(rowmax > .nonFuncThreshold)
      # funcRegIndex <- which(rowSums(RegWeightMat) > nonFuncThreshold)
      RegulatoryCoordinateExtractorOutput <- RegulatoryCoordinateExtractorOutput[funcRegIndex, ]
      RegWeightMat <- RegWeightMat[funcRegIndex, ]
    }
  }
    
  # if removeNonFunctionalRegs option is on, modify RegulatoryCoordinateExtractorOutput, RegWeightMat and RowSideCol
  if(removeNonFunctionalRegs){
    result <- cbind(RegulatoryCoordinateExtractorOutput, RegWeightMat)
  }
  return(result)
}


####################################################################################################################################################
####################################################################################################################################################
#example
#a <- RegulatoryCoordinateExtractor(gene.Names = rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2))
#aa <- EnhWeightReader("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_17/Cluster_1/Outputs")
aaaa <- sort(Exp17.GEMSTAT.output.sigmoid_3_softmax_5$averageSSE, index.return = T, decreasing = F)$ix
aaa <- RegElemWeightCoordinate(RegulatoryCoordinateExtractorOutput = a,
                               EnhWeightReaderOutput = aa,
                               modelIndex = aaaa[1:19], removeNonFunctionalRegs = F, nonFuncThreshold = 0.5)

####################################################################################################################################################
####################################################################################################################################################
ParameterHeatMap <- function(parameterMat,
                             exportplot = T,
                             filename = "parameterHeatMap.png",
                             .Rowv = T,
                             .Colv = F,
                             .dendrogram = "row",
                             logTransform = T,
                             .cellnote = format(round(log10(parameterMat), 2), nsmall = 2),
                             .RowSideColors = character(0),
                             col_breaks = numeric(0),
                             col_break_ch_nu = 3,
                             rownames_care=T,
                             .rowsep=integer(0),
                             .colsep=integer(0)
                             #nu_colors = 299
                             ){
  # this function draws a heatmap of the log of parameters of multiple models
  # parameterMat is the matrix of parameter values, where each row is a model and each column is a parameter
  # .Rowv : determines if and how the row dendrogram should be reordered.
  # .dendrogram : character string indicating whether to draw 'none', 'row', 'column' or 'both' dendrograms.
  # .RowSideColors : character vector with length equal to number of models. specifies the color in side bar for each model
  # logTransform : boolean, if True the input matrix will be log transformed before visualization
  # col_breaks : the numbers to move to the next color in the map
  # col_break_ch_nu : number of different base colors to be used
  # rownames_care : if True, cares about the rownames of parameterMat and uses them other wise doesn't care
  # .rowsep :  vector of integers indicating which rows should be separated from the preceding columns or rows by a narrow space of color sepcolor.
  # .colsep :  vector of integers indicating which columns should be separated from the preceding columns or rows by a narrow space of color sepcolor.
  
  if(rownames_care){
    rownames(parameterMat) <- unlist(lapply(strsplit(rownames(parameterMat), split = "\\."),"[[", 2))
  }
  library(RColorBrewer)
  library(gplots)
  if(length(.RowSideColors) == 0){
    .RowSideColors = rep("white", nrow(parameterMat))
  }
  
  if (length(col_breaks)==0){
    col_breaks = c(seq(-4, -2, length=100),     # for red
                   seq(-1.99, 0, length=100),           # for yellow
                   seq(0.01, parRange[2], length=100)) 
  }
  
  all_my_cols <- c("white", "blue", "cyan", "green", "yellow", "red", "black")
  my_cols <- all_my_cols[(4 - floor((col_break_ch_nu - 1)/2)):(4 + ceiling((col_break_ch_nu - 1)/2))]
  my_palette <- colorRampPalette(my_cols)(n = (length(col_breaks)-1))
  #export or not?
  if(exportplot){
    png(filename,    # create PNG for the heat map        
        width = 8*300,        # 5 x 300 pixels
        height = 8*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
  }
  if(logTransform){
    parameterMat <- log10(parameterMat)
  }
  #plot the heatmap using gplot heatmap.2
  parRange <- range(parameterMat)
  

            # for green
  heatmap.2(parameterMat,
            #cellnote = .cellnote,  # same data set for cell labels
            main = "parameter log scale", # heat map title
            notecol="black",      # change font color of cell labels to black
            # density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(9,7),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            breaks=col_breaks,    # enable color transition at specified limits
            dendrogram=.dendrogram,     # only draw a col dendrogram
            RowSideColors = .RowSideColors, # grouping row-variables into different categories
            Rowv = .Rowv,
            Colv = .Colv,
            scale = "none",
            rowsep = .rowsep,
            colsep = .colsep
  ) 
  if(exportplot){
    dev.off() 
  }
}
####################################################################################################################################################
####################################################################################################################################################
#example
par(mfrow = c(1,1), mar = c(6,4,4,6))
aaaa <- sort(Exp17.GEMSTAT.output.sigmoid_3_softmax_5$averageSSE, index.return = T, decreasing = F)$ix
ParameterHeatMap(parameterMat = Exp17.GEMSTAT.param.perf.sigmoid_3_softmax_5$parameters[aaaa[1:19], ])
####################################################################################################################################################
####################################################################################################################################################
PerformanceHeattMap <- function(GEMSTAT_output_Reader_ensemble.Output,
                                geneNames = rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2),
                                modelIndex = c(1:length(GEMSTAT_output_Reader_ensemble.Output$Allresults[[1]])),
                                geneIndex = c(1:length(GEMSTAT_output_Reader_ensemble.Output$Allresults[[1]][[1]]$SSE)),
                                .Colv = T, .Rowv = F, .dendrogram = "col",exportplot = T, filename = "PerformanceHeatMap.png",
                                .RowSideColors = character(0), .ColSideColors = character(0)){
  # creates a heatmap of performance of all chosen models
  # GEMSTAT_output_Reader_ensemble.Output is the output of GEMSTAT_output_Reader_ensemble function
  # .Rowv : determines if and how the row dendrogram should be reordered.
  # .Colv : determines if and how the column dendrogram should be reordered.
  # .dendrogram : character string indicating whether to draw 'none', 'row', 'column' or 'both' dendrograms.
  # .RowSideColors is a charachter vector with length equal to number of genes indicating the color of the sidebar for each gene
  # .ColSideColors is a charachter vector with length equal to number of models indicating the color of the sidebar for each model
  #number of genes
  if(length(.RowSideColors) == 0){
    .RowSideColors = rep("white", length(geneIndex))
  }
  if(length(.ColSideColors) == 0){
    .ColSideColors = rep("white", length(modelIndex))
  }
  ngenes <- length(geneIndex)
  #create the SSE matrix and name rows and columns
  SSEmat <- matrix(1L, nrow = ngenes, ncol = length(modelIndex))
  rownames(SSEmat) <- geneNames[geneIndex]
  colnames(SSEmat) <- unlist(lapply(strsplit(names(GEMSTAT_output_Reader_ensemble.Output$Allresults[[1]][modelIndex]), split = "\\."),"[[", 2))
  for(model in 1:length(modelIndex)){
    SSEmat[, model] <- GEMSTAT_output_Reader_ensemble.Output$Allresults[[1]][[modelIndex[model]]]$SSE[geneIndex]
  }
  #colors for heatmap
  my_palette <- colorRampPalette(c("green", "yellow", "red"))(n = 299)
  #export or not?
  if(exportplot){
    png(filename,    # create PNG for the heat map        
        width = 8*300,        # 5 x 300 pixels
        height = 8*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
  }
  #plot the heatmap using gplot heatmap.2
  
  col_breaks = c(seq(0.0, 0.149, length=100),     # for green
                 seq(0.150, 0.199, length=100),   # for yellow
                 seq(0.2, 0.6, length=100))         # for red
  heatmap.2(SSEmat,
            #cellnote = format(round(SSEmat, 2), nsmall = 2),  # same data set for cell labels
            main = "performance", # heat map title
            notecol="black",      # change font color of cell labels to black
            # density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(9,7),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            breaks=col_breaks,    # enable color transition at specified limits
            dendrogram=.dendrogram,     # only draw a .dendrogram dendrogram
            RowSideColors = .RowSideColors,# grouping row-variables(genes) into different categories
            ColSideColors = .ColSideColors,# grouping column-variables(models) into different categories
            Rowv = .Rowv,
            Colv = .Colv,
            scale = "none"
            ) 
  if(exportplot){
    dev.off() 
  }
}
####################################################################################################################################################
####################################################################################################################################################
#example
aad2 <- sort(Exp17.GEMSTAT.output.sigmoid_no_softmax_5$averageSSE, decreasing = F, index.return = T)$ix 
PerformanceHeattMap(Exp17.GEMSTAT.output.sigmoid_no_softmax_5, modelIndex = c(1:200))
####################################################################################################################################################
####################################################################################################################################################

CompareModelPlot <- function(GEMSTAT_output_Reader_ensemble.Output,
                             LogFileReader.Output,
                             EnhWeightReader.Output,
                             modelIndex,
                             ..col_vec = col_vector,
                             .geneNames=rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2),
                             ..exportplot = T,
                             filename = paste("CompareModels",modelIndex[1], modelIndex[2], sep = "_"),
                             plotDir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Ensemble/EXP17",
                             .setmaxto1RestZero = F){
  # This function creates a plot that compare two models in the following way:
  # compare log10(parameters) (heatmap)
  # compute the sum of squared differences in Contribution of reg elements for each gene between the two models: plot a barplot for all genes
  # plot heatmap of enhancer contribution for the genes that are different in their contribution (see above)
  # plot Histogram of difference in performance for all genes between the two models
  # plot Expression profile of some of the genes that were predicted better than a threshold in at least one of the models
  # inputs: Outputs of three functions: GEMSTAT_output_Reader_ensemble, LogFileReader , EnhWeightReader
  # modelIndex : is an integer vector of length 2 containing the index of the models to compare
  # ..col_vec is a vector of color names
  # .geneNames : charachtar vector containing the names of the genes
  # .setmaxto1RestZero :  if True for each gene sets the weight of the enhancer with max weight to 1 and rest of its enhancers to zero (used depending on the GEMSTAT obj func)
  
  # Creates total of six plots
  #                      heatmap of parameters,
  #                      hist of difference in performance for all genes between the two models
  #                      barplot of difference in enhancer contrib,
  #                      heatmap of enhancer contribution for the genes that are different in their contribution,
  #                      heat map of model performance
  #                       plot Expression profile of some of the genes 
  
  

  # layout(matrix(c(1,2,3,4,4,4,5,5,5), 3, 3, byrow = F), 
  #        widths=c(2,1,1), heights = c(0.8,1,1))
  prev_dir <- getwd()
  setwd(plotDir)
  dir.create(filename)
  setwd(filename)
  # first heatmap of parameters
  ParameterHeatMap(LogFileReader.Output$parameters[modelIndex,], exportplot = T, filename = "parameterHeatMap.png", .Rowv = F, .dendrogram = "none")
  # second hist of difference in performance
  perfDif <- GEMSTAT_output_Reader_ensemble.Output$Allresults[[1]][[modelIndex[1]]]$SSE - GEMSTAT_output_Reader_ensemble.Output$Allresults[[1]][[modelIndex[2]]]$SSE
  if(..exportplot){
    png("PerfDif.png",    # create PNG for the heat map
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
  }
  hist(perfDif,breaks = 60, col = 3, main = "SSE (model1) - SSE (model2)", xlab = "")
  if(..exportplot){
    dev.off()
  }
  #Third hist difference in enhancer contribution
  EnhWeights <- EnhWeightReader.Output[modelIndex]
  #function to compute SSE between enhancer contribution lists
  SSEListCompute <- function(List1, List2){
    result <- numeric(length(List1))
    for (i in 1:length(result)){
      result[i] <- sqrt(sum((List1[[i]] - List2[[i]])^2)/length(List1[[i]]))
    }
    return(result)
  }
  EnhWeightsSSE <- SSEListCompute(EnhWeights[[1]], EnhWeights[[2]])
  if(..exportplot){
    png("EnhWeightsSSE.png",    # create PNG for the heat map
        width = 8*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
  }
  barplot(EnhWeightsSSE, col = 3, main = "Contributing Reg difference")
  if(..exportplot){
    dev.off()
  }
  #Forth plot: Enh Contribution heatmap
  RegElemWeightHeatMap(.EnhWeightReaderOutput = EnhWeightReader.Output,
                       .modelIndex = modelIndex,
                       .removeNonFunctionalRegs = F,
                       .exportplot = T,
                       .col_vec = ..col_vec,
                       .dendrogram = "none",
                       .Colv = F,
                       setmaxto1RestZero = .setmaxto1RestZero)
  #Fifth plot: performance heatmap
  PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output,
                      modelIndex = modelIndex,
                      geneNames = .geneNames,
                      .Colv = F,.Rowv = F, .dendrogram = "none", exportplot = T)
  #Sixth plot: performance heatmap
  # get the prediction of the models for all genes
  cur_Prediction1 <- GEMSTAT_output_Reader_ensemble.Output$Allresults[[1]][[modelIndex[1]]]$prediction
  cur_Prediction2 <- GEMSTAT_output_Reader_ensemble.Output$Allresults[[1]][[modelIndex[2]]]$prediction
  
  # get the true value of expression for all genes
  cur_real <- GEMSTAT_output_Reader_ensemble.Output$Allresults[[1]][[1]]$TrueValue
  rownames(cur_real) <- .geneNames
  # get the SSE of each gene
  cur_SSE1 <- GEMSTAT_output_Reader_ensemble.Output$Allresults[[1]][[modelIndex[1]]]$SSE
  cur_SSE2 <- GEMSTAT_output_Reader_ensemble.Output$Allresults[[1]][[modelIndex[2]]]$SSE
  #filter genes that have high SSE in both models
  FilterGenes <- (! rowSums(cbind((cur_SSE1 > 0.2), (cur_SSE2 > 0.2))) %in% 2)
  cur_Prediction1 <- cur_Prediction1[FilterGenes, ]
  cur_Prediction2 <- cur_Prediction2[FilterGenes, ]
  cur_real <- cur_real[FilterGenes, ]
  cur_SSE1 <- cur_SSE1[FilterGenes]
  cur_SSE2 <- cur_SSE2[FilterGenes]
  if(..exportplot){
    png("predictionvsReal.png",    # create PNG for the heat map
        width = 11*300,        # 5 x 300 pixels
        height = 11*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
  }
  par(mfrow = c(12, 11), mar = c(0.05, 0.05, 0.05, 0.05))
  for(gene in 1:nrow(cur_Prediction1)){#plot for each gene
    plot(cur_real[gene,],
         ylim = c(0,1), col = 1,
         xaxt = "n", yaxt = "n", main ="", ylab = "", xlab = "",
         type = "l", lwd = 2)
    
    lines(cur_Prediction1[gene,], ylim = c(0,1), col = 2)
    lines(cur_Prediction2[gene,], ylim = c(0,1), col = 3)
    aa1 <- cur_SSE1[gene]
    aa2 <- cur_SSE2[gene]
    text((ncol(cur_Prediction1) - 2), 0.9 , format(round(aa1, 2), nsmall = 2), cex = 0.5) # write the size of the SSE
    text((ncol(cur_Prediction2) - 2), 0.8 , format(round(aa2, 2), nsmall = 2), cex = 0.5) # write the size of the SSE
    
    text(2, 0.1 , rownames(cur_real)[gene], cex = 0.5) # write the index of the gene

  }
  if(..exportplot){
    dev.off()
  }

  setwd(prev_dir)
}

####################################################################################################################################################
####################################################################################################################################################
#example
CompareModelPlot(GEMSTAT_output_Reader_ensemble.Output = Exp17.GEMSTAT.output.sigmoid_no_softmax_5,
                 LogFileReader.Output = Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5,
                 EnhWeightReader.Output = Exp17.GEMSTAT.enhWeight.sigmoid_no_softmax_5,
                 modelIndex = c(1,2),
                 ..col_vec = col_vector,
                 .geneNames=rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2),
                 ..exportplot = T,
                 plotDir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Ensemble/EXP17")
####################################################################################################################################################
####################################################################################################################################################
PredictedvsRealCorRMSD <- function(GEMSTAT_output_Reader_ensemble.output,
                                   modelIndex = c(1:length(GEMSTAT_output_Reader_ensemble.output$Allresults[[1]])),
                                   geneIndex = c(1:nrow(GEMSTAT_output_Reader_ensemble.output$Allresults[[1]][[1]]$prediction)),
                                   corrMethod = "both"){
  # this function computes the (pearson and spearman) correlation of real vs predicted for
  # each gene in each model, given the output of GEMSTAT_output_Reader_ensemble function
  # corrMethod can be either "pearson", "spearman", or "both"
  # modelIndex is integer vector containing the index of models to be investigated
  # geneIndex is integer vector containing the index of genes to be investigated
  # output is a list with three entries:
   # first entry is the matrix of pearson corr : nrow : number of genes, ncol : number of models
   # second entry is the matrix of spearman corr : dimensions same as first entry
   # Third entry is the matrix of RMSD which is just extracted from input: dimensions same as first entry
  
  nmodels <- length(modelIndex)
  ngenes <- length(geneIndex)
  corMatPearson <- matrix(nrow = ngenes, ncol = nmodels) #matrix that holds pearson correlation
  rownames(corMatPearson) <- geneIndex
  colnames(corMatPearson) <- modelIndex
  corMatSpearman <- matrix(nrow = ngenes, ncol = nmodels) #matrix that holds spearman correlation
  rownames(corMatSpearman) <- geneIndex
  colnames(corMatSpearman) <- modelIndex
  RMSDMat <- matrix(nrow = ngenes, ncol = nmodels)  #matrix that holds RMSD values
  rownames(RMSDMat) <- geneIndex
  colnames(RMSDMat) <- modelIndex
  for(model in 1:nmodels){
    if (corrMethod != "spearman"){
      corMatPearson[, model] <- myCorTest(GEMSTAT_output_Reader_ensemble.output$Allresults[[1]][[modelIndex[model]]]$prediction[geneIndex, ],
                                          GEMSTAT_output_Reader_ensemble.output$Allresults[[1]][[modelIndex[model]]]$TrueValue[geneIndex, ],
                                          lenOu = ngenes,
                                          CorMethod = "pearson")[[2]]
    }
    if (corrMethod != "pearson"){
      corMatSpearman[, model] <- myCorTest(GEMSTAT_output_Reader_ensemble.output$Allresults[[1]][[modelIndex[model]]]$prediction[geneIndex, ],
                                          GEMSTAT_output_Reader_ensemble.output$Allresults[[1]][[modelIndex[model]]]$TrueValue[geneIndex, ],
                                          lenOu = ngenes,
                                          CorMethod = "spearman")[[2]]
    }
    RMSDMat[, model] <- GEMSTAT_output_Reader_ensemble.output$Allresults[[1]][[modelIndex[model]]]$SSE[geneIndex]
  }
  if (corrMethod == "both"){
    result = list(pearson = corMatPearson, spearman = corMatSpearman, RMSD = RMSDMat)
  }else if(corrMethod == "pearson"){
    result = list(pearson = corMatPearson, RMSD = RMSDMat)
  }else{
    result = list(spearman = corMatSpearman, RMSD = RMSDMat)
  }
  return(result)
}

####################################################################################################################################################
####################################################################################################################################################
#example
aa <- PredictedvsRealCorRMSD(GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.All.Filtered$output,corrMethod = "both")
par(mfrow = c(3,1), mar = c(1,2,3,1))
boxplot.matrix(aa$pearson[,aad2], main = "Pearson Correlation")
boxplot.matrix(aa$spearman[,aad2], main = "Spearman Correlation")
boxplot.matrix(aa$RMSD[,aad2], main = "RMSD")

####################################################################################################################################################
####################################################################################################################################################
PlotPercentGenes <- function(..GEMSTAT_output_Reader_ensemble.output, .LogFileReader.output, .EnhWeightReader.output,
                             .geneIndex = c(1:nrow(..GEMSTAT_output_Reader_ensemble.output$Allresults[[1]][[1]]$prediction)),
                             .modelIndex = c(1:length(..GEMSTAT_output_Reader_ensemble.output$Allresults[[1]])),
                             RMSDrange = c(0,0.15,0.2,1),
                             CorRange = c(1, 0.8, 0.6, -1),
                             .corrMethod = "pearson",
                             Num.Cluster = 4,
                             sortBasedOn = 1,
                             .CorRMSDList=numeric(0),
                             exportplot = F){
  # This function plots the percent of genes that have RMSD or Correlation between some ranges
  # inputs:
  # .GEMSTAT_output_Reader_ensemble.output : output of .GEMSTAT_output_Reader_ensemble.output
  # .LogFileReader.output is the ouput of LogFileReader function
  # .EnhWeightReader.output is the ouput of EnhWeightReader function
  # .geneIndex : index of the genes that we care about
  # .modelIndex : index of the models that we care about
  # RMSDrange defines the ranges that barplots show for RMSD
  # CorRange defines the ranges that barplots show for Cor
  # .corrMethod is either "pearson" or "spearman"
  # Num.Cluster is the number of clusters
  # sortBasedOn is the range which you want to sort with respect to
  # .CorRMSDList is the output of PredictedvsRealCorRMSD function for the same inputs: be aware that if this is given then .modelIndex, .geneIndex, .corrMethod should have been considred in it.
  
  # compute the correlations and store RMSD in one matrix
  if(length(.CorRMSDList) == 0){
    CorRMSDList <- PredictedvsRealCorRMSD(GEMSTAT_output_Reader_ensemble.output = ..GEMSTAT_output_Reader_ensemble.output,
                                          modelIndex = .modelIndex,
                                          geneIndex = .geneIndex,
                                          corrMethod = .corrMethod )
  }else{
    CorRMSDList <- .CorRMSDList
  }

  #create matrices to store number of genes in each range 
  CorRangeMat <- matrix(nrow = 3, ncol = length(.modelIndex))
  # rownames(CorRange) <- c(paste(CorRange[1]),
  #                         paste(),
  #                         paste())
  RMSDRangeMat <- matrix(nrow = 3, ncol = length(.modelIndex))
  rownames(RMSDRangeMat) <- character(nrow(RMSDRangeMat))
  colnames(RMSDRangeMat) <- .modelIndex
  rownames(CorRangeMat) <- character(nrow(CorRangeMat))
  colnames(CorRangeMat) <- .modelIndex
  for (rng in 1:(length(RMSDrange) - 1)){
    aa1 <- CorRMSDList[[2]] >= RMSDrange[rng]
    aa2 <- CorRMSDList[[2]] < RMSDrange[rng + 1]
    aa = aa1 * aa2
    RMSDRangeMat[rng, ] <- colSums(aa)/length(.geneIndex)
    rownames(RMSDRangeMat)[rng] <- paste(RMSDrange[rng], "> RMSD >",RMSDrange[rng+1])
    
    cc1 <- CorRMSDList[[1]] <= CorRange[rng]
    cc2 <- CorRMSDList[[1]] > CorRange[rng + 1]
    cc = cc1 * cc2
    CorRangeMat[rng, ] <- colSums(cc)/length(.geneIndex)
    rownames(CorRangeMat)[rng] <- paste(CorRange[rng], "> Cor >",CorRange[rng+1])

  } # end of counting genes explained by each model
  # sort based on the number of genes based explained
  RMSDRangeMat.sort.index <- sort(RMSDRangeMat[sortBasedOn,], decreasing = T, index.return = T)$ix
  CorRangeMat.sort.index  <- sort(CorRangeMat[sortBasedOn,] , decreasing = T, index.return = T)$ix
  RMSDRangeMat.Sorted <- RMSDRangeMat[, RMSDRangeMat.sort.index]
  CorRangeMat.Sorted <- CorRangeMat[, CorRangeMat.sort.index]
  cat("CorRangeMat.Sorted is \n")
  print(CorRangeMat.Sorted[,1:10])
  ######################## cluster based on parameter values ########################
  paramCluster <- kmeans(log10(.LogFileReader.output$parameters[.modelIndex, ]), Num.Cluster)
  paramClusterIndexRMSDSorted <- integer(0)
  paramClusterIndexCorSorted <- integer(0)
  paramClusterBorder <- integer(0)
  for(cl in 1:Num.Cluster){
    clInd <- which(paramCluster$cluster %in% cl)
    #sort based in best explaining models RMSD
    RMSDsortedInd <- sort(RMSDRangeMat[sortBasedOn, clInd], decreasing = T, index.return = T)$ix 
    paramClusterIndexRMSDSorted <- c(paramClusterIndexRMSDSorted, clInd[RMSDsortedInd])
    #sort based in best explaining models Cor
    CorsortedInd <- sort(CorRangeMat[sortBasedOn, clInd], decreasing = T, index.return = T)$ix 
    paramClusterIndexCorSorted <- c(paramClusterIndexCorSorted, clInd[CorsortedInd])
    #border assignment
    paramClusterBorder <- c(paramClusterBorder,sum(paramCluster$size) - sum(paramCluster$size[(Num.Cluster - cl + 1):Num.Cluster]))
  }
  
  ######################## cluster based on enh contributions ########################
  #create a matrix containing all regulatory elements contribution: each row a reg elem, each column a model
  RegElNumber <- sum(unlist(lapply(.EnhWeightReader.output[[1]][.geneIndex], length)))
  RegWeightMat <- matrix(nrow = RegElNumber, ncol = length(.modelIndex))
  for(model in 1:length(.modelIndex)){
    RegWeightMat[, model] <- unlist(.EnhWeightReader.output[[.modelIndex[model]]][.geneIndex])
  }
  colnames(RegWeightMat) <- as.character(.modelIndex)
  #cluster based on enhancer contribution
  EnhContrCluster <- kmeans(t(RegWeightMat), Num.Cluster)
  EnhContrClusterIndexRMSDSorted <- integer(0)
  EnhContrClusterIndexCorSorted <- integer(0)
  #EnhContrClusterBorder <- 0
  for(cl in 1:Num.Cluster){
    clInd <- which(EnhContrCluster$cluster %in% cl)
    #sort based in best explaining models RMSD
    RMSDsortedInd <- sort(RMSDRangeMat[sortBasedOn, clInd], decreasing = T, index.return = T)$ix 
    EnhContrClusterIndexRMSDSorted <- c(EnhContrClusterIndexRMSDSorted, clInd[RMSDsortedInd])
    #sort based in best explaining models Cor
    CorsortedInd <- sort(CorRangeMat[sortBasedOn, clInd], decreasing = T, index.return = T)$ix 
    EnhContrClusterIndexCorSorted <- c(EnhContrClusterIndexCorSorted, clInd[CorsortedInd])
  }
  # EnhContrClusterBorder <- integer(Num.Cluster + 1)
  # EnhContrClusterBorder[1] <- 0
  # for(i in 1:(Num.Cluster)){
  #   EnhContrClusterBorder[i+1] <- EnhContrClusterBorder[i] + EnhContrCluster$size[i]
  # }
  ################################## plot ############################################
  prev.Par <- par() #save the previous par setting
  if(exportplot){
    png("PercentExplained_RMSD.png",    # create PNG for the heat map
        width = 11*300,        # 5 x 300 pixels
        height = 11*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
  }
  par(mfrow = c(3, 1), mar = c(2,3,4,1))
  #plotting RMSD
  barplot(RMSDRangeMat.Sorted, xlab = "", main = "RMSD all", border=NA, las=2)
  abline(h = seq(0.1, 1, 0.1) , col = 2, lwd = 0.2, lty = 2)
  barplot(RMSDRangeMat[, paramClusterIndexRMSDSorted], xlab = "", main = "parameter clustered", border=NA, las=2)
  #abline(v = paramClusterBorder , col = 2)
  abline(h = seq(0.1, 1, 0.1) , col = 2, lwd = 0.2, lty = 2)
  barplot(RMSDRangeMat[, EnhContrClusterIndexRMSDSorted], xlab = "", main = "Enh contrib clustered", border=NA, las=2)
  #abline(v = EnhContrClusterBorder , col = 2)
  abline(h = seq(0.1, 1, 0.1) , col = 2, lwd = 0.2, lty = 2)
  if(exportplot){
    dev.off()
  }
  
  #plotting Cor
  if(exportplot){
    png("PercentExplained_Cor.png",    # create PNG for the heat map
        width = 11*300,        # 5 x 300 pixels
        height = 11*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
  }
  par(mfrow = c(3, 1), mar = c(2,3,4,1))
  barplot(CorRangeMat.Sorted, xlab = "", main = "Correlation all", border=NA, las=2)
  abline(h = seq(0.1, 1, 0.1) , col = 2, lwd = 0.2, lty = 2)
  barplot(CorRangeMat[, paramClusterIndexCorSorted], xlab = "", main = "parameter clustered", border=NA, las=2)
  abline(h = seq(0.1, 1, 0.1) , col = 2, lwd = 0.2, lty = 2)
  #abline(v = paramClusterBorder , col = 2)
  barplot(CorRangeMat[, EnhContrClusterIndexCorSorted], xlab = "", main = "Enh contrib clustered", border=NA, las=2)
  #abline(v = EnhContrClusterBorder, col = 2)
  abline(h = seq(0.1, 1, 0.1) , col = 2, lwd = 0.2, lty = 2)
  
  if(exportplot){
    dev.off()
  }
  par(mfrow = prev.Par$mfrow, mar = prev.Par$mar)
  Clusterings <- list(paramBased = paramCluster, enhBased = EnhContrCluster)
  PercentExp <- list(RMSD = RMSDRangeMat, Cor = CorRangeMat)
  result = list(PercentExplained = PercentExp, Clustering = Clusterings)
  return(result)
}

####################################################################################################################################################
####################################################################################################################################################
#example
aaPR <- PredictedvsRealCorRMSD(GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.All.Filtered$output,
                               corrMethod = "pearson" )

aaa <-      PlotPercentGenes(..GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.All.Filtered$output,
                             .LogFileReader.output=Exp17.GEMSTAT.All.Filtered$param,
                             .EnhWeightReader.output=Exp17.GEMSTAT.All.Filtered$enh,
                             RMSDrange = c(0,0.15,0.2,1),
                             CorRange = c(1, 0.8, 0.6, -1),
                             .corrMethod = "pearson",
                             Num.Cluster = 4,
                             sortBasedOn = 1,
                             .CorRMSDList = aaPR)
####################################################################################################################################################
####################################################################################################################################################
PurturbationInputCreator <- function(modelParameters,
                                     TFnames=c("ESR1", "RARA"),
                                     .coopertingTFIndex= rbind(c(1,1),c(1,2)),
                                     modelIndex,
                                     experiment.nu,
                                     .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                     .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                                     .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/"){
  
  #Gets the model parameters (a numeric vector of length 7) in the order output by logfilereader function: 1.ESR1_Binding  2.ESR1_txp 3.RARA_binding  4.RARA_txp   5.basal 6.ESR1_ESR1 7.RARA_ESR1
  # creates parameter files for purturbation experiments
  # creates jobs for purturbation experiments
  # copies everything to hal
  
  # create a purturbation parameter folder if it doesn't exist
  # create a purturbation output folder if it doesn't exist
  prev_dir <- getwd()
  setwd(dir = .dir)
  dir.create(paste0(.dir,"/Experiment_", experiment.nu, "/Cluster_1/Inputs/PurturbParameters"))
  dir.create(paste0(.dir,"/Experiment_", experiment.nu, "/Cluster_1/Outputs/PurturbOutputs"))
  
  setwd(paste0(.dir,"/Experiment_", experiment.nu, "/Cluster_1/Inputs/PurturbParameters"))
  #create 4 purturbed versions of the parameter vector: #1 not purturbed #2 set ESR1_txp to 1, #3 set RARA_txp to 1 ,#4 set ESR1_ESR1 to 1 ,#5 set RARA_ESR1 to 1
  purturbedModels <- matrix(nrow = 5, ncol = length(modelParameters))
  purturbedModels[1, ] <- modelParameters
  parChan <- c(2, 4, 6, 7)
  for (i in 2:5){
    purturbedModels[i, ] <- modelParameters
    purturbedModels[i, parChan[i-1]] <- 1
  }
  for(samp in 1:nrow(purturbedModels)){
    curSample <- purturbedModels[samp, ]
    for(TFnu in 1:length(TFnames)){
      #write TF parameters
      cat(paste(TFnames[TFnu], curSample[(TFnu)*2 -1],
                curSample[(TFnu)*2], sep = "\t"),
          file = paste0("par_purt_",modelIndex ,"_", samp,".txt"),
          sep = "\n", append = T)
    }
    #write basal
    cat(paste0("basal_transcription = ", curSample[(TFnu)*2 + 1]),
        file = paste0("par_purt_", modelIndex, "_", samp,".txt"), sep = "\n", append = T)
    #write coop
    for(cptf in 1:nrow(.coopertingTFIndex)){
      cat(paste(TFnames[.coopertingTFIndex[cptf,1]],
                TFnames[.coopertingTFIndex[cptf,2]],
                curSample[length(TFnames)*2 + 1 + cptf], sep = "\t"),
          file = paste0("par_purt_", modelIndex, "_", samp,".txt"), sep = "\n", append = T)
    }
  }#wrote the parameterFiles
  ########  ########  ########  ########  ########  ########  ########
  setwd(paste0(.dir,"/Experiment_", experiment.nu))
  #create the job file
  CurJob <- readLines(paste0("Experiment_",experiment.nu,"job"))[1]
  for(purt in 1:nrow(purturbedModels)){
    tmpJob <- gsub(pattern = "src_TrySoftMax", replacement = "src_Max", x = CurJob)
    tmpJob <- gsub(pattern =  paste0("/Outputs/Experiment_",experiment.nu,"_Cluster_1\\.outpar","\\S*","\\.txt\\.outEns"),
                   replacement = paste0("/Outputs/PurturbOutputs/Experiment_",experiment.nu,"_",modelIndex,"_",purt,".outEns"),
                   x = tmpJob)
    tmpJob <- gsub(pattern = paste0("Parameters/par","\\S*","\\.txt"),
                   replacement = paste0("PurturbParameters/par_purt_",modelIndex,"_",purt,".txt"),
                   x = tmpJob)
    tmpJob <- gsub(pattern = "-na 10", replacement = "-na 0", x = tmpJob)
    tmpJob <- gsub(pattern = paste0("/Outputs/Experiment_",experiment.nu,"_Cluster_1.logpar","\\S*",".txt.logEns"),
                   replacement = paste0("/Outputs/PurturbOutputs/Experiment_",experiment.nu,"_",modelIndex,"_",purt,".logEns"),
                   x = tmpJob)
    cat(tmpJob,
        file = paste0("Experiment_", "Purturb_",experiment.nu,"_",modelIndex,"job"), sep = "\n", append = T)
  }
  ########  ########  ########  ########  ########  ########  ########
  # create a .submit from the written hal job
  sys_com <- paste0("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl"," ",paste0("Experiment_", "Purturb_",experiment.nu,"_",modelIndex,"job "),
                    5," tmpjob_", experiment.nu, "_", modelIndex,
                    " >job_purturb_exp_", experiment.nu, "_", modelIndex, ".submit")
  cat(c("#!/bin/bash", "\n"),file=paste0("hal_sub_creator_purturb_", "exp_", experiment.nu, "_", modelIndex), sep="")
  cat(sys_com, file=paste0("hal_sub_creator_purturb_", "exp_", experiment.nu, "_", modelIndex), sep="", append = T)
  
  sys_com <- paste0("chmod +x ", paste0(" hal_sub_creator_purturb_", "exp_", experiment.nu, "_", modelIndex))
  system(sys_com)
  #####################################################
  #copy to hal
  #create the bash file
  cat(c("#!/bin/bash", "\n"),file=paste0("PurtJob_", "exp_", experiment.nu,"_",modelIndex), sep="")
  # mount thr directory you want to remote:
  cat(c(paste0("sshfs"," tabebor2@veda.cs.illinois.edu:" ,.sharedmountsdir,
               " ~/remote -oauto_cache,reconnect,defer_permissions,noappledouble,negative_vncache,volname=MySSHFSMount"), "\n"),
      file=paste0("PurtJob_", "exp_", experiment.nu,"_",modelIndex), sep="", append=T)
  # cd to ~/remote directory
  cat(c("cd ~/remote", "\n"),
      file=paste0("PurtJob_", "exp_", experiment.nu,"_",modelIndex), sep="", append=T)
  source_dir <- paste0(.dir,"/Experiment_", experiment.nu, "/Cluster_1/Inputs/PurturbParameters")
  cat(c(paste("cp -r ", source_dir, " ", "Experiment_", experiment.nu, "/Cluster_1/Inputs/" , sep=""), "\n"),
      file=paste0("PurtJob_", "exp_", experiment.nu,"_",modelIndex),
      sep="", append=T)
  source_dir <- paste0(.dir,"/Experiment_", experiment.nu, "/Cluster_1/Outputs/PurturbOutputs")
  cat(c(paste("cp -r ", source_dir, " Experiment_", experiment.nu,"/Cluster_1/Outputs/" , sep=""), "\n"),
      file=paste0("PurtJob_", "exp_", experiment.nu,"_",modelIndex),
      sep="", append=T)
  source_dir <- paste0(.dir,"/Experiment_", experiment.nu)
  cat(c(paste("cp ", source_dir,
              paste0("/hal_sub_creator_purturb_", "exp_", experiment.nu, "_", modelIndex),
              " Experiment_", experiment.nu,"/" , sep=""), "\n"),
      file=paste0("PurtJob_", "exp_", experiment.nu,"_",modelIndex),
      sep="", append=T)
  cat(c(paste("cp ", source_dir, paste0("/Experiment_", "Purturb_",experiment.nu,"_",modelIndex,"job")," Experiment_", experiment.nu,"/" , sep=""), "\n"),
      file=paste0("PurtJob_", "exp_", experiment.nu,"_",modelIndex),
      sep="", append=T)
  cat(c("cd ~", "\n"),
      file=paste0("PurtJob_", "exp_", experiment.nu,"_",modelIndex), sep="", append=T)
  cat(c("umount", "remote/", "\n"),
      file=paste0("PurtJob_", "exp_", experiment.nu,"_",modelIndex), sep=" ", append=T)
  #make the created bash file executable
  sys_com <- paste0("chmod +x ", .dir,"/Experiment_", experiment.nu,
                    paste0("/PurtJob_", "exp_", experiment.nu,"_",modelIndex))
  system(sys_com)
  print("Copying to veda ...")
  sys_com <- paste0(.dir,"/Experiment_", experiment.nu,
                    paste0("/PurtJob_", "exp_", experiment.nu,"_",modelIndex))
  system(sys_com)
  
  setwd(prev_dir)
}
###############################################################################################################
###############################################################################################################
#example
PurturbationInputCreator(modelParameters = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[14,],
                         TFnames=c("ESR1", "RARA"),
                         .coopertingTFIndex= rbind(c(1,1),c(1,2)),
                         modelIndex=14,
                         experiment.nu = 17,
                         .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                         .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                         .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/")
####################################################################################################################################################
####################################################################################################################################################
Purturbation_output_Reader <- function(.exp.nu,
                                       modelIndex,
                                       .root_dir="tabebor2@veda.cs.illinois.edu:/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble",
                                       .dest_dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble"
                                       ,copyOutput = T){
  # 1) .exp.nu : number of the experiment
  # 2) modelIndex = Index of the model: this is actually not important since it extracts output of all models, but names the bash file after this model
  # 3) .root_dir : the root directory in veda where Experiment folders are located
  # 4) .dest_dir : the destination directory on my laptop where experiment folders are located.
  
  prev_dir <- getwd()
  setwd(paste0(.dest_dir, "/Experiment_", .exp.nu))
  if(copyOutput){
    cat(c("#!/bin/bash", "\n"),file=paste0("bash_directory_output_exp_", .exp.nu,"_",modelIndex), sep="")
    # mount the directory you want to remote:
    cat(c(paste("sshfs", .root_dir,
                "~/remote -oauto_cache,reconnect,defer_permissions,noappledouble,negative_vncache,volname=MySSHFSMount"), "\n"),
        file=paste0("bash_directory_output_exp_", .exp.nu,"_",modelIndex),
        sep="", append=T)
    # cd to ~/remote directory
    cat(c("cd ~/remote", "\n"),
        file=paste0("bash_directory_output_exp_", .exp.nu,"_",modelIndex),
        sep="", append=T)
    cat(c(paste("cp ","Experiment_", .exp.nu,"/Cluster_",1,"/Outputs/PurturbOutputs/","* ",
                .dest_dir, "/Experiment_", .exp.nu,"/Cluster_",1,"/Outputs/PurturbOutputs/", sep=""), "\n"),
        file=paste0("bash_directory_output_exp_", .exp.nu,"_",modelIndex),
        sep="", append=T)
    cat(c("cd ~", "\n"),
        file=paste0("bash_directory_output_exp_", .exp.nu,"_",modelIndex),
        sep="", append=T)
    cat(c("umount", "remote/", "\n"),
        file=paste0("bash_directory_output_exp_", .exp.nu,"_",modelIndex),
        sep=" ", append=T)
    # make the created bash file executable
    sys_com <- paste0("chmod +x ",paste0("bash_directory_output_exp_", .exp.nu,"_",modelIndex))
    system(sys_com)
    print("Copying from veda ...")
    sys_com <- paste0("./", paste0("bash_directory_output_exp_", .exp.nu,"_",modelIndex))
    system(sys_com)
  }
  #reading the copied output
  AllResults <- list()
  tmpFile <- list.files(paste0(.dest_dir, "/Experiment_",
                               .exp.nu,"/Cluster_",1,"/Outputs/PurturbOutputs/"),pattern = "*.outEns")
  par_number <- length(tmpFile)
  SSEResults <- matrix(nrow =par_number , ncol = 1)
  rownames(SSEResults) <- tmpFile
  for(clus in 1:1){
    AllResults[[clus]] <- list()
    print(paste0("reading output for cluster_",clus,"..."))
    output_file_names <- list.files(paste0(.dest_dir, "/Experiment_",
                                           .exp.nu,"/Cluster_",clus,"/Outputs/PurturbOutputs/"),pattern = "*.outEns")
    OutPut_clus <- list() #list contaning the outputs of this cluster
    for(output_file in 1:length(output_file_names)){
      OutPut_clus[[output_file]] <- read.table(paste0(.dest_dir, "/Experiment_",
                                                      .exp.nu,"/Cluster_",clus,"/Outputs/PurturbOutputs/", output_file_names[output_file]),
                                               header = T, sep = "\t")
      cur_Prediction <- OutPut_clus[[output_file]][seq(2, nrow(OutPut_clus[[output_file]]), 2), 2:ncol(OutPut_clus[[output_file]])]
      cur_real <- OutPut_clus[[output_file]][seq(1, nrow(OutPut_clus[[output_file]]), 2), 2:ncol(OutPut_clus[[output_file]])]
      cur_Prediction <- do.call(cbind, cur_Prediction)
      cur_real <- do.call(cbind, cur_real)
      rownames(cur_Prediction) <- unique(OutPut_clus[[output_file]][, 1])
      rownames(cur_real) <- rownames(cur_Prediction)
      cur_SSE <- sqrt(rowSums((cur_real - cur_Prediction)^2)/ncol(cur_real))
      names(cur_SSE) <- rownames(cur_real) 
      AllResults[[clus]][[output_file]] <- list(prediction=cur_Prediction, TrueValue=cur_real, SSE = cur_SSE)
      SSEResults[output_file, clus] <- mean(cur_SSE)
    }
    names(AllResults[[clus]]) <- output_file_names
  }#end of loop over clusters
  setwd(prev_dir)
  return(list(Allresults = AllResults, averageSSE = SSEResults))
}
#######################################################################################################################
#######################################################################################################################
#example
aaa <- Purturbation_output_Reader (.exp.nu = 17,
                                   modelIndex = 14,
                                   .root_dir="tabebor2@veda.cs.illinois.edu:/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble",
                                   .dest_dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble"
                                   ,copyOutput = T)
####################################################################################################################################################
####################################################################################################################################################

####################################################################################################################################################
####################################################################################################################################################
#example

####################################################################################################################################################
####################################################################################################################################################
#next function starts here

####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
###################################      ENSEMBLE EXPERIMENTS      #################################################################################   
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################

######################   EXP 11   ###############################################################
#creating ranges for parameter sampling
aa <- matrix(0L, nrow = 6, ncol = 3)
aa[1,] <- c(-2, 4, 1)
aa[2,] <- c(-2, 4, 1)
aa[3,] <- c(-2, 4, 1)
aa[4,] <- c(-2, 4, 1)
aa[5,] <- c(-2, 3, 1)
aa[6,] <- c(-2, 3, 1)
#creating parameter samples, setting up gemstat
GEMSTAT_input_constructor_Ensemble(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed))),
                                   gene_expression_Mat = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed, 
                                   .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                   .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                   .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                   .promoter.sequence=Promoter.Sequence.Char,
                                   .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                   .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                                   .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                   experiment.nu = 11,
                                   .motifs=TF.motifs.TimeSeries.count,
                                   .TF_expression=TF.Expression.TimeSeries,
                                   .max.enh.length = numeric(0),
                                   ..o="DIRECT", coopertingTFIndex=rbind(c(1,1),c(1,2))
                                   , ..oo="SSE", ..ct=13, ..sigma=character(0)
                                   ,.par.RangeMat = aa, .nu.sample= 20)
#reading gemstat results
Exp11.GEMSTAT.output <- GEMSTAT_output_Reader_ensemble(.exp.nu = 11, .number.Of.Clusters = 1)
#plotting best results
par(mfrow = c(1,1), mar = c(4,4,4,4))
hist(Exp11.GEMSTAT.output$averageSSE[,1], breaks = 100)
sort(Exp11.GEMSTAT.output$averageSSE[,1],decreasing = F ,index.return = T)$ix
plot_ensemble_output(Exp11.GEMSTAT.output, ClusterNumber = 1, paramnumber = 26)


aa <- matrix(0L, nrow = length(Exp11.GEMSTAT.output$Allresults[[1]]), ncol = nrow(Exp11.GEMSTAT.output$Allresults[[1]][[1]]$prediction))
for (i in 1:length(Exp11.GEMSTAT.output$Allresults[[1]])){
  for (j in 1:nrow(Exp11.GEMSTAT.output$Allresults[[1]][[i]]$prediction)){
    aa[i,j] <- var(Exp11.GEMSTAT.output$Allresults[[1]][[i]]$prediction[j,4:10])
  }
}
aaa <- apply(aa, 1, max)
sort(aaa,decreasing = T, index.return = T)$ix
par(mfrow = c(1,1), mar = c(4,4,4,4))
boxplot.matrix(aa, use.cols = F)
###look to find how many of these genes are form the ER-RAR associated group and how many from ER alone group, plot separately
rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed)
sum(rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm) %in% Genes.Associated.REMAP.ER.Entrez.TS.DIff.NoRAR.Sampled)
sum(rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm) %in% Genes.Associated.RAR.ER.Entrez.TS.DIff)
sum(rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed) %in% Genes.Associated.REMAP.ER.Entrez.TS.DIff.NoRAR.Sampled)
sum(rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed) %in% Genes.Associated.RAR.ER.Entrez.TS.DIff)

sum(rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.EarlyPeak) %in% Genes.Associated.REMAP.ER.Entrez.TS.DIff.NoRAR.Sampled)
sum(rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.EarlyPeak) %in% Genes.Associated.RAR.ER.Entrez.TS.DIff)
aa <- GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed[rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed) %in% Genes.Associated.RAR.ER.Entrez.TS.DIff, ]
aaa <- GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed[rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed) %in% Genes.Associated.REMAP.ER.Entrez.TS.DIff.NoRAR.Sampled, ]
par(mfrow = c(9,9), mar = c(0.1, 0.1, 0.1, 0.1))

for(i in 1:nrow(aa)){
  plot(aa[i,], main = "", ylim = c(0,1), xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l")
}
par(mfrow = c(7,6), mar = c(0.1, 0.1, 0.1, 0.1))

for(i in 1:nrow(aaa)){
  plot(aaa[i,], main = "", ylim = c(0,1), xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l")
}



######################   EXP 12   ###############################################################
#Same as experiment 11 but filtering out the genes that have significant increased expression after the third time point

aa <- apply(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed, 1, which.max)
GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.EarlyPeak <- GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed[aa < 5, ]


#creating ranges for parameter sampling
aa <- matrix(0L, nrow = 6, ncol = 3)
aa[1,] <- c(-2, 4, 1)
aa[2,] <- c(-2, 4, 1)
aa[3,] <- c(-2, 4, 1)
aa[4,] <- c(-2, 4, 1)
aa[5,] <- c(-2, 3, 1)
aa[6,] <- c(-2, 3, 1)
#creating parameter samples, setting up gemstat
GEMSTAT_input_constructor_Ensemble(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.EarlyPeak))),
                                   gene_expression_Mat = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.EarlyPeak, 
                                   .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                   .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                   .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                   .promoter.sequence=Promoter.Sequence.Char,
                                   .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                   .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                                   .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                   experiment.nu = 12,
                                   .motifs=TF.motifs.TimeSeries.count,
                                   .TF_expression=TF.Expression.TimeSeries,
                                   .max.enh.length = numeric(0),
                                   ..o="DIRECT", coopertingTFIndex=rbind(c(1,1),c(1,2))
                                   , ..oo="SSE", ..ct=13, ..sigma=character(0)
                                   ,.par.RangeMat = aa, .nu.sample= 20)
#reading gemstat results
Exp12.GEMSTAT.output <- GEMSTAT_output_Reader_ensemble(.exp.nu = 12, .number.Of.Clusters = 1)

#plotting best results
par(mfrow = c(1,1), mar = c(4,4,4,4))
hist(Exp12.GEMSTAT.output$averageSSE[,1], breaks = 100)
sort(Exp12.GEMSTAT.output$averageSSE[,1],decreasing = F ,index.return = T)$ix
plot_ensemble_output(Exp12.GEMSTAT.output, ClusterNumber = 1, paramnumber = 221)


aa <- matrix(0L, nrow = length(Exp12.GEMSTAT.output$Allresults[[1]]), ncol = nrow(Exp12.GEMSTAT.output$Allresults[[1]][[1]]$prediction))
for (i in 1:length(Exp12.GEMSTAT.output$Allresults[[1]])){
  for (j in 1:nrow(Exp12.GEMSTAT.output$Allresults[[1]][[i]]$prediction)){
    aa[i,j] <- var(Exp12.GEMSTAT.output$Allresults[[1]][[i]]$prediction[j,4:10])
  }
}
aaa <- apply(aa, 1, median)
sort(aaa,decreasing = T, index.return = T)$ix
par(mfrow = c(1,1), mar = c(4,4,4,4))
boxplot.matrix(aa, use.cols = F)

names(Exp12.GEMSTAT.output$Allresults[[1]])[221]


######################   EXP 13   ###############################################################

####Try to find a new subset of ERonly genes in which early repressed and late enhanced patterns are not present
Genes.Associated.REMAP.ER.Entrez.TS.DIff.NoRAR.Sampled2 <- sample(x = setdiff(Genes.Associated.REMAP.ER.Entrez.TS.DIff,
                                                                             Genes.Associated.RAR.ER.Entrez.TS.DIff), size = 535, replace = F)
Genes.Associated.REMAP.ER.Entrez.TS.DIff.NoRAR.Sampled2 <- setdiff(Genes.Associated.REMAP.ER.Entrez.TS.DIff.NoRAR.Sampled2, aaChosen)
a <- match(Genes.Associated.REMAP.ER.Entrez.TS.DIff.NoRAR.Sampled2, rownames(GSE78167.RNAseq.Avg))
aa <- GSE78167.RNAseq.Avg.Norm[a, ]
par(mfrow = c(10,10), mar = c(0.1, 0.1, 0.1, 0.1))
for(i in 1:nrow(aa)){
  plot(aa[i,], main = "", ylim = c(0,1), xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l")
  text(2, 0.9 , i, cex = 0.5) # write the size of the SSE
  
}

a1 <- c(9,11,17,25,42,43,48,76,92,97,189,228,235,236,239,249,253,292,320,331,334,361,373,377,381)

aaChosen <- unique(c(aaChosen, Genes.Associated.REMAP.ER.Entrez.TS.DIff.NoRAR.Sampled2[a1]))
a <- match(aaChosen, rownames(GSE78167.RNAseq.Avg.Norm))
aa2 <- GSE78167.RNAseq.Avg.Norm[a, ]



aa <- GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed[rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed) %in% Genes.Associated.RAR.ER.Entrez.TS.DIff, ]
aa <- aa[-c(6,10,14,16,22,24,30,40,47,53,54,62,66),]


GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered <- rbind(aa, aa2)

#creating ranges for parameter sampling
aa <- matrix(0L, nrow = 6, ncol = 3)
aa[1,] <- c(-2, 4, 1)
aa[2,] <- c(-2, 4, 1)
aa[3,] <- c(-2, 4, 1)
aa[4,] <- c(-2, 4, 1)
aa[5,] <- c(-2, 3, 1)
aa[6,] <- c(-2, 3, 1)
#creating parameter samples, setting up gemstat
GEMSTAT_input_constructor_Ensemble(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered))),
                                   gene_expression_Mat = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered, 
                                   .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                   .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                   .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                   .promoter.sequence=Promoter.Sequence.Char,
                                   .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                   .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                                   .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                   experiment.nu = 13,
                                   .motifs=TF.motifs.TimeSeries.count,
                                   .TF_expression=TF.Expression.TimeSeries,
                                   .max.enh.length = numeric(0),
                                   ..o="DIRECT", coopertingTFIndex=rbind(c(1,1),c(1,2))
                                   , ..oo="SSE", ..ct=13, ..sigma=character(0)
                                   ,.par.RangeMat = aa, .nu.sample= 20)

#reading gemstat results
Exp13.GEMSTAT.output <- GEMSTAT_output_Reader_ensemble(.exp.nu = 13, .number.Of.Clusters = 1)

#plotting best results
par(mfrow = c(1,1), mar = c(4,4,4,4))
hist(Exp13.GEMSTAT.output$averageSSE[,1], breaks = 100)
sort(Exp13.GEMSTAT.output$averageSSE[,1],decreasing = F ,index.return = T)$ix
plot_ensemble_output(Exp13.GEMSTAT.output, ClusterNumber = 1, paramnumber = 111)


aa <- matrix(0L, nrow = length(Exp13.GEMSTAT.output$Allresults[[1]]), ncol = nrow(Exp13.GEMSTAT.output$Allresults[[1]][[1]]$prediction))
for (i in 1:length(Exp13.GEMSTAT.output$Allresults[[1]])){
  for (j in 1:nrow(Exp13.GEMSTAT.output$Allresults[[1]][[i]]$prediction)){
    aa[i,j] <- var(Exp13.GEMSTAT.output$Allresults[[1]][[i]]$prediction[j,4:10])
  }
}
aaa <- apply(aa, 1, median)
sort(aaa,decreasing = T, index.return = T)$ix
par(mfrow = c(1,1), mar = c(4,4,4,4))
boxplot.matrix(aa, use.cols = F)

names(Exp13.GEMSTAT.output$Allresults[[1]])[196]


######################   EXP 14   ###############################################################
#repeating experiment 13, but storing log files to work with parameters
#creating ranges for parameter sampling
aa <- matrix(0L, nrow = 6, ncol = 3)
aa[1,] <- c(-2, 4, 1)
aa[2,] <- c(-2, 4, 1)
aa[3,] <- c(-2, 4, 1)
aa[4,] <- c(-2, 4, 1)
aa[5,] <- c(-2, 3, 1)
aa[6,] <- c(-2, 3, 1)
#creating parameter samples, setting up gemstat
GEMSTAT_input_constructor_Ensemble(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered))),
                                   gene_expression_Mat = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered, 
                                   .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                   .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                   .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                   .promoter.sequence=Promoter.Sequence.Char,
                                   .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                   .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                                   .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                   experiment.nu = 14,
                                   .motifs=TF.motifs.TimeSeries.count,
                                   .TF_expression=TF.Expression.TimeSeries,
                                   .max.enh.length = numeric(0),
                                   ..o="DIRECT", coopertingTFIndex=rbind(c(1,1),c(1,2))
                                   , ..oo="SSE", ..ct=13, ..sigma=character(0)
                                   ,.par.RangeMat = aa, .nu.sample= 20)

#reading gemstat results
Exp14.GEMSTAT.output <- GEMSTAT_output_Reader_ensemble(.exp.nu = 14, .number.Of.Clusters = 1)
plot_ensemble_output(Exp14.GEMSTAT.output, ClusterNumber = 1, paramnumber = 134)




Exp14.GEMSTAT.param.perf <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_14/Cluster_1/Outputs")
par(mforw = c(1,1), mar = c(4,4,4,4))
hist(Exp14.GEMSTAT.param.perf$performance, breaks = 100)
aa2 <- kmeans(Exp14.GEMSTAT.param.perf$parameters, 4)
aa22 <- aa2$cluster
aa <- dist(Exp14.GEMSTAT.param.perf$parameters)
aaa <- cmdscale(aa,eig=TRUE, k=2)

# plot solution 
a1 <- aaa$points[,1]
a2 <- aaa$points[,2]
plot(a1, a2, xlab="Coordinate1", ylab="Coordinate2", 
     main="Metric MDS", col = aa22 + 1)


par(mforw = c(1,1), mar = c(7,3,4,4))
for(i in 1:4){
  boxplot.matrix(log10(Exp14.GEMSTAT.param.perf$parameters)[aa2$cluster == i,], ylim = c(-5,4),main=paste0("cl ",i," , mean perf: ", format (round(mean(Exp14.GEMSTAT.param.perf$performance[aa2$cluster == i]), 4), nsmall = 4)), las =2)
  text(7, -4 , sum(aa2$cluster == i), cex = 0.5) # write the size of the SSE
}


#sort genes based on the variance of performance in different models
plot_ensemble_output(Exp14.GEMSTAT.output, ClusterNumber=1, OneGeneManyModel = T, GeneNumber = 0, VarRank = 1)
##looking at models with ER-ER > 10 and ER-RAR > 10
aa <- intersect(which(log10(Exp14.GEMSTAT.param.perf$parameters[,6]) > 1),
                which(log10(Exp14.GEMSTAT.param.perf$parameters[,7]) > 1))
##looking at models with ER-ER > 10 and ER-RAR < 1
aaa <- intersect(which(log10(Exp14.GEMSTAT.param.perf$parameters[,6]) > 1) ,
                 which(log10(Exp14.GEMSTAT.param.perf$parameters[,7]) < 0)) 
##looking at models with ER-ER > 10 and ER-RAR > 10 and ER txp > 1 and RAR txp > 1
aaaa <- intersect(aa, intersect(which(log10(Exp14.GEMSTAT.param.perf$parameters[,2]) > 0),
                                which(log10(Exp14.GEMSTAT.param.perf$parameters[,4]) > 0)))
##looking at models with ER-ER > 10 and ER-RAR > 10 and ER txp > 1 and RAR txp > 1 : None
aaaaa <- intersect(aa, intersect(which(log10(Exp14.GEMSTAT.param.perf$parameters[,4]) > 0),
                                 which(log10(Exp14.GEMSTAT.param.perf$parameters[,4]) < 0)))
##looking at models with ER-ER > 1 and ER-RAR > 1
aa <- intersect(which(log10(Exp14.GEMSTAT.param.perf$parameters[,6]) > 0),
                which(log10(Exp14.GEMSTAT.param.perf$parameters[,7]) > 0))
##looking at models with ER-ER > 1 and ER-RAR < 1
aaa <- intersect(which(log10(Exp14.GEMSTAT.param.perf$parameters[,6]) > 0),
                 which(log10(Exp14.GEMSTAT.param.perf$parameters[,7]) < 0)) 
##looking at models with ER-ER > 1 and ER-RAR > 1 and ER txp > 1 and RAR txp > 1 
aaaa <- intersect(aa, intersect(which(log10(Exp14.GEMSTAT.param.perf$parameters[,2]) > 0),
                                which(log10(Exp14.GEMSTAT.param.perf$parameters[,4]) > 0)))
##looking at models with ER-ER > 1 and ER-RAR > 1 and ER txp > 1 and RAR txp < 1 : only 1: in which it essentially doesn't use RAR
aaaa <- intersect(aa, intersect(which(log10(Exp14.GEMSTAT.param.perf$parameters[,2]) > 0),
                                which(log10(Exp14.GEMSTAT.param.perf$parameters[,4]) < 0)))
##looking at models with ER-ER > 10 and ER-RAR > 1
aa <- intersect(which(log10(Exp14.GEMSTAT.param.perf$parameters[,6]) > 1),
                which(log10(Exp14.GEMSTAT.param.perf$parameters[,7]) > 0))
##looking at models with ER-ER > 10 and ER-RAR < 1
aaa <- intersect(which(log10(Exp14.GEMSTAT.param.perf$parameters[,6]) > 1),
                 which(log10(Exp14.GEMSTAT.param.perf$parameters[,7]) < 0)) 

aa4 <- kmeans(Exp14.GEMSTAT.param.perf$parameters, 4)
##looking at models in cluster1 vs cluster 2
aa <- which(aa4$cluster %in% 1)
aaa <- which(aa4$cluster %in% 2)
plot_ensemble_output(Exp14.GEMSTAT.output,
                     ClusterNumber=1,
                     OneGeneManyModel = T, Rank = 2,
                     groupvsgroup = T,
                     group1 = aa, group2 = aaa)

##looking at models in cluster1 vs cluster 3
aa <- which(aa4$cluster %in% 1)
aaa <- which(aa4$cluster %in% 3)
plot_ensemble_output(Exp14.GEMSTAT.output,
                     ClusterNumber=1,
                     OneGeneManyModel = T, Rank = 2,
                     groupvsgroup = T,
                     group1 = aa, group2 = aaa)
##looking at models in cluster1 vs cluster 4
##looking at models in cluster2 vs cluster 3
##looking at models in cluster2 vs cluster 4
##looking at models in cluster3 vs cluster 4






#plot_ensemble_output(Exp14.GEMSTAT.output, ClusterNumber=1, OneGeneManyModel = T, GeneNumber = 0, VarRank = 1, varAmong = c(aa,aaa), plotAmong = c(aa,aaa))

plot_ensemble_output(Exp14.GEMSTAT.output,
                                 ClusterNumber=1, paramnumber=0,
                                 OneGeneManyModel = T, GeneNumber = 0, Rank = 2, varAmong = integer(0), plotAmong = integer(0),
                                 groupvsgroup = T, group1 = aa, group2 = aaa)
# 
# 
# aa <- matrix(nrow = length(Exp14.GEMSTAT.output$Allresults[[1]][[1]]$SSE), ncol = length(Exp14.GEMSTAT.output$Allresults[[1]]))
# for (i in 1:length(Exp14.GEMSTAT.output$Allresults[[1]])){
#   aa[,i] <- Exp14.GEMSTAT.output$Allresults[[1]][[i]]$SSE
# }
# 
# aaa <- apply(aa,1,sd)
# sort(aaa, decreasing = T, index.return = T)$ix
# 


######################   EXP 15   ###############################################################
#removing genes not explanable using current TF profiles
GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2 <- GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered[-c(88), ]

aa <- matrix(0L, nrow = 6, ncol = 3)
aa[1,] <- c(-3, 2, 1) #ESR1_Binding
aa[2,] <- c(0, 4, 1) #ESR1_txp
aa[3,] <- c(-3, 2, 1) #RARA_binding
aa[4,] <- c(-2, 4, 1) #RARA_txp
aa[5,] <- c(-2, 4, 1) #ESR1_ESR1
aa[6,] <- c(-2, 4, 1) #RARA_ESR1
#creating parameter samples, setting up gemstat
GEMSTAT_input_constructor_Ensemble(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2))),
                                   gene_expression_Mat = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2, 
                                   .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                   .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                   .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                   .promoter.sequence=Promoter.Sequence.Char,
                                   .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                   .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                                   .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                   experiment.nu = 15,
                                   .motifs=TF.motifs.TimeSeries.count,
                                   .TF_expression=TF.Expression.TimeSeries,
                                   .max.enh.length = numeric(0),
                                   ..o="DIRECT", coopertingTFIndex=rbind(c(1,1),c(1,2))
                                   , ..oo="SSE", ..ct=11, ..sigma=character(0)
                                   ,.par.RangeMat = aa, .nu.sample= 50)


#reading gemstat results
Exp15.GEMSTAT.output <- GEMSTAT_output_Reader_ensemble(.exp.nu = 15, .number.Of.Clusters = 1, copyOutput = F)
aa <- sort(Exp15.GEMSTAT.output$averageSSE, decreasing = F, index.return = T)$ix 
plot_ensemble_output(Exp15.GEMSTAT.output, ClusterNumber = 1, paramnumber = 1119)
plot_ensemble_output(Exp15.GEMSTAT.output, ClusterNumber = 1, OneGeneManyModel = T, Rank = 1)
plot_ensemble_output(Exp15.GEMSTAT.output, ClusterNumber = 1, OneGeneManyModel = T, GeneNumber = 86)

Exp15.GEMSTAT.param.perf <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_15/Cluster_1/Outputs", removeDups = F)


par(mforw = c(1,1), mar = c(4,4,4,4))
hist(Exp15.GEMSTAT.param.perf$performance, breaks = 100)
aa2 <- kmeans(Exp15.GEMSTAT.param.perf$parameters, 5)
aa22 <- aa2$cluster
aa <- dist(Exp15.GEMSTAT.param.perf$parameters)
aaa <- cmdscale(aa,eig=TRUE, k=2)

# plot solution 
a1 <- aaa$points[,1]
a2 <- aaa$points[,2]
plot(a1, a2, xlab="Coordinate1", ylab="Coordinate2", 
     main="Metric MDS", col = aa22 + 1)


par(mforw = c(1,1), mar = c(7,3,4,4))
for(i in 1:5){
  boxplot.matrix(log10(Exp15.GEMSTAT.param.perf$parameters)[aa2$cluster %in% i,], ylim = c(-5,4),main=paste0("cl ",i," , mean perf: ", format (round(mean(Exp15.GEMSTAT.param.perf$performance[aa2$cluster == i]), 4), nsmall = 4)), las =2)
  text(7, -4 , sum(aa2$cluster == i), cex = 0.5) # write the size of the SSE
}
sum(log10(Exp15.GEMSTAT.param.perf$parameters[,6]) > 0)
hist(Exp15.GEMSTAT.param.perf$performance[log10(Exp15.GEMSTAT.param.perf$parameters[,6]) > 0])

aa <- which(log10(Exp15.GEMSTAT.param.perf$parameters[,6]) > 1)
length(aa)
plot_ensemble_output(Exp15.GEMSTAT.output, ClusterNumber = 1, OneGeneManyModel = T, Rank = 3,varAmong = aa,plotAmong = aa)


intersect(which((log10(Exp15.GEMSTAT.param.perf$parameters[,6]) > 1)), which(Exp15.GEMSTAT.param.perf$performance < 0.186))
which((log10(Exp15.GEMSTAT.param.perf$parameters[,6]) > 1) && which(Exp15.GEMSTAT.param.perf$performance < 0.19))

Exp15.GEMSTAT.param.perf$performance[log10(Exp15.GEMSTAT.param.perf$parameters[,6]) > 1 ]

aa <- intersect(which((log10(Exp15.GEMSTAT.param.perf$parameters[,6]) > 1)), which(Exp15.GEMSTAT.param.perf$performance < 0.186))
aa1 <- intersect(aa, which(log10(Exp15.GEMSTAT.param.perf$parameters[,7]) < 0))


boxplot.matrix(log(Exp15.GEMSTAT.param.perf$parameters)[aa,])
##############################################################
######################   EXP 16   ###############################################################
#removing genes not explanable using current TF profiles
GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2 <- GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered[-c(88), ]

aa <- matrix(0L, nrow = 6, ncol = 3)
aa[1,] <- c(-4, 2, 1) #ESR1_Binding
aa[2,] <- c(0, 4, 1) #ESR1_txp
aa[3,] <- c(-4, 2, 1) #RARA_binding
aa[4,] <- c(-2, 4, 1) #RARA_txp
aa[5,] <- c(0, 4, 1) #ESR1_ESR1
aa[6,] <- c(-2, 4, 1) #RARA_ESR1
#creating parameter samples, setting up gemstat
GEMSTAT_input_constructor_Ensemble(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2))),
                                   gene_expression_Mat = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2, 
                                   .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                   .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                   .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                   .promoter.sequence=Promoter.Sequence.Char,
                                   .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                   .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                                   .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                   experiment.nu = 16,
                                   .motifs=TF.motifs.TimeSeries.count,
                                   .TF_expression=TF.Expression.TimeSeries,
                                   .max.enh.length = numeric(0),
                                   ..o="DIRECT", coopertingTFIndex=rbind(c(1,1),c(1,2))
                                   , ..oo="SSE", ..ct=12, ..sigma=character(0)
                                   ,.par.RangeMat = aa, .nu.sample= 100)


#reading gemstat results
Exp16.GEMSTAT.output <- GEMSTAT_output_Reader_ensemble(.exp.nu = 16, .number.Of.Clusters = 1, copyOutput = F)
Exp16.GEMSTAT.param.perf <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_16/Cluster_1/Outputs", removeDups =F)



aad <- sort(Exp16.GEMSTAT.param.perf$performance, decreasing = F, index.return = T)$ix 
plot_ensemble_output(Exp16.GEMSTAT.output, ClusterNumber = 1, paramnumber = 2708)
plot_ensemble_output(Exp16.GEMSTAT.output, ClusterNumber = 1, OneGeneManyModel = T, Rank = 1)
plot_ensemble_output(Exp16.GEMSTAT.output, ClusterNumber = 1, OneGeneManyModel = T, GeneNumber = 86)

par(mforw = c(1,1), mar = c(4,4,4,4))
hist(Exp16.GEMSTAT.param.perf$performance, breaks = 100)
aa2 <- kmeans(Exp16.GEMSTAT.param.perf$parameters, 5)
aa22 <- aa2$cluster
aa <- dist(Exp16.GEMSTAT.param.perf$parameters)
aaa <- cmdscale(aa,eig=TRUE, k=2)

# plot solution 
a1 <- aaa$points[,1]
a2 <- aaa$points[,2]
aapch <- rep(0, length(a2))
aapch[aa[1:300]] <- rep(4, 300)
plot(a1, a2, xlab="Coordinate1", ylab="Coordinate2", 
     main="Metric MDS", col = aa22 + 1, pch = aapch)
points(a1[aa[1:300]], a2[aa[1:300]], cex = 3)

par(mforw = c(1,1), mar = c(7,3,4,4))
for(i in 1:5){
  boxplot.matrix(log10(Exp16.GEMSTAT.param.perf$parameters)[aa2$cluster %in% i,],ylim = c(-11,9),
                 main=paste0("cl ",i," , mean perf: ",
                             format (round(mean(Exp16.GEMSTAT.param.perf$performance[aa2$cluster == i]), 4), nsmall = 4)), las =2)
  text(7, -4 , sum(aa2$cluster == i), cex = 0.5) # write the size of the SSE
}


sum(log(Exp16.GEMSTAT.param.perf$parameters)[,6] > 0)

aa <- which(log(Exp16.GEMSTAT.param.perf$parameters)[,6] > 1)
aaa <- which(log(Exp16.GEMSTAT.param.perf$parameters)[,7] >0 )

length(intersect(aa, aaa))
hist(Exp16.GEMSTAT.param.perf$performance[intersect(aa, aaa)], breaks = 100)
hist(Exp16.GEMSTAT.param.perf$performance[setdiff(aa, aaa)], breaks = 100)

plot_ensemble_output(Exp16.GEMSTAT.output, ClusterNumber = 1, paramnumber = 1119)

boxplot.matrix(log(Exp16.GEMSTAT.param.perf$parameters[aa[1:500], ]))

#######
#filter for models which ER-ER interaction is greater than 1
aaERER <- which(log10(Exp16.GEMSTAT.param.perf$parameters)[,6] > 0)
aaERRAR <- which(log10(Exp16.GEMSTAT.param.perf$parameters)[,7] > 0)
aaERERRAR <- intersect(aaERER, aaERRAR)
aaERERNoRAR <- setdiff(aaERER, aaERRAR)
plot_ensemble_output(Exp16.GEMSTAT.output, ClusterNumber = 1, paramnumber = aaERERRAR[1])
plot_ensemble_output(Exp16.GEMSTAT.output, ClusterNumber = 1, paramnumber = aaERERNoRAR[1])

aa <- PerformanceFilter(ensembleOutput = Exp16.GEMSTAT.output ,cluster = 1,
                        modelIndex=numeric(0), geneIndex = c(16, 13),
                        perfThresh = 0.18)


#cluster the real expression values:

aa <- clustrator(ExpressionMat = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2, nfolds = 4,exportPlot = F)
Genes.Clustered.Based.On.Expression <- aa
aaa <- c(sample(aa$ClusterIndex[[1]], size = length(aa$ClusterIndex[[2]]), replace = F), aa$ClusterIndex[[2]])
aa2 <- PerformanceFilter(ensembleOutput = Exp16.GEMSTAT.output ,cluster = 1,
                        modelIndex=numeric(0), geneIndex = aa$ClusterIndex[[1]], sortModels = T, strictRemove = F)
plot_ensemble_output(Exp16.GEMSTAT.output, ClusterNumber = 1, paramnumber = aa2[1])


length(intersect(aa2[1:500], aaERERRAR))
par(mar = c(7,4,4,4))
boxplot.matrix(log10(Exp16.GEMSTAT.param.perf$parameters[aa2[1:500],]), las = 2, main = "logscale")


par(mfrow = c(4,1), mar = c(2,3,1,1))
aa2 <- PerformanceFilter(ensembleOutput = Exp16.GEMSTAT.output ,cluster = 1,
                         modelIndex=numeric(0), geneIndex = aa$ClusterIndex[[1]], sortModels = T, strictRemove = F)
boxplot.matrix(log10(Exp16.GEMSTAT.param.perf$parameters[aa2[1:100],]), las = 2, ylim = c(-5,4), xaxt = "n")
abline(h = 0, col = 2)



aa2 <- PerformanceFilter(ensembleOutput = Exp16.GEMSTAT.output ,cluster = 1,
                         modelIndex=numeric(0), geneIndex = aa$ClusterIndex[[2]], sortModels = T, strictRemove = F)
boxplot.matrix(log10(Exp16.GEMSTAT.param.perf$parameters[aa2[1:100],]), las = 2, ylim = c(-5,4),  xaxt = "n")
abline(h = 0, col = 2)

aa2 <- PerformanceFilter(ensembleOutput = Exp16.GEMSTAT.output ,cluster = 1,
                         modelIndex=numeric(0), geneIndex = aa$ClusterIndex[[3]], sortModels = T, strictRemove = F)
boxplot.matrix(log10(Exp16.GEMSTAT.param.perf$parameters[aa2[1:100],]), las = 2, ylim = c(-5,4),  xaxt = "n")
abline(h = 0, col = 2)

par(mar = c(4,3,1,1))
aa2 <- PerformanceFilter(ensembleOutput = Exp16.GEMSTAT.output ,cluster = 1,
                         modelIndex=numeric(0), geneIndex = aa$ClusterIndex[[4]], sortModels = T, strictRemove = F)
boxplot.matrix(log10(Exp16.GEMSTAT.param.perf$parameters[aa2[1:100],]), las = 1, ylim = c(-5,4))
abline(h = 0, col = 2)




sum(log10(Exp16.GEMSTAT.param.perf$parameters[,7]) > 0)
sum(log10(Exp16.GEMSTAT.param.perf$parameters[,6]) > 0)
sum(log10(Exp16.GEMSTAT.param.perf$parameters[,7]) < 0)

boxplot(log10(Exp16.GEMSTAT.param.perf$parameters[,7]) )


aaa <- c(sample(aa$ClusterIndex[[3]], size = length(aa$ClusterIndex[[2]]), replace = F), aa$ClusterIndex[[2]])
aa2 <- PerformanceFilter(ensembleOutput = Exp16.GEMSTAT.output ,cluster = 1,
                         modelIndex=numeric(0), geneIndex = c(aa$ClusterIndex[[1]], aa$ClusterIndex[[4]]), sortModels = T, strictRemove = F)

par(mfrow = c(1,1), mar = c(4,4,4,4))
boxplot.matrix(log10(Exp16.GEMSTAT.param.perf$parameters), las = 1, ylim = c(-5,4))
abline(h = 0, col = 2)



Exp16.GEMSTAT.param.perf$parameters[aa2[1], ]
###################
# plot the number of well modeled genees in each subgroup of models
par(mfrow = c(4,1), mar = c(2,3,1,1))
ff <- integer(0)
for(cl in 1:4){
  aa2 <- PerformanceFilter(ensembleOutput = Exp16.GEMSTAT.output ,cluster = 1,
                           modelIndex=numeric(0), geneIndex = Genes.Clustered.Based.On.Expression$ClusterIndex[[cl]], sortModels = T, strictRemove = F)
  aaa <- matrix(nrow = nrow(aa2), ncol = nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2))
  for(i in 1:nrow(aa2)){
    aaa[i,] <- Exp16.GEMSTAT.output$Allresults[[1]][[aa2[i,1]]]$SSE
  }
  aaam <- apply(aaa[1:200, ], 2, median)
  #sum(aaam < 0.15)
  boxplot.matrix(aaa[1:200, ], las = 2, ylim = c(0,0.4), xaxt = "n", outline=FALSE)
  abline(h = 0.15, col = 2)
  text(60, 0.02 , paste0(sum(aaam < 0.15), " out of ", ncol(aaa) ," genes: ", "median SSE less than 0.15"), cex = 1) # write the size of the SSE
  ff <- union(ff, (which(aaam < 0.15)))
}


aaa <- c(sample(Genes.Clustered.Based.On.Expression$ClusterIndex[[3]], size = length(Genes.Clustered.Based.On.Expression$ClusterIndex[[2]]), replace = F), Genes.Clustered.Based.On.Expression$ClusterIndex[[2]])
aa2 <- PerformanceFilter(ensembleOutput = Exp16.GEMSTAT.output ,cluster = 1,
                         modelIndex=numeric(0), geneIndex = aaa, sortModels = T, strictRemove = F)
aaa <- matrix(nrow = nrow(aa2), ncol = nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2))
for(i in 1:nrow(aa2)){
  aaa[i,] <- Exp16.GEMSTAT.output$Allresults[[1]][[aa2[i,1]]]$SSE
}
aaam <- apply(aaa[1:500, ], 2, median)
par(mfrow = c(1,1), mar = c(4,4,4,4))

boxplot.matrix(aaa[1:500, ], las = 2, ylim = c(0,0.4), xaxt = "n", outline=FALSE)

abline(h = 0.15, col = 2)
text(60, 0.02 , paste0(sum(aaam < 0.15), " out of ", ncol(aaa) ," genes: ", "median SSE less than 0.15"), cex = 1) # write the size of the SSE

par(mfrow = c(7,6), mar = c(0.1,0.1,0.1,0.1))
for(i in 1:length(intersect(which(aaam < 0.15), ff))){
  plot(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2[i,],
       ylim = c(0,1), col = 2,
       xaxt = "n", yaxt = "n", main ="", ylab = "", xlab = "",
       type = "l", lwd = 2)
}


######################   EXP 16_using_sigmoid 1/(1 + exp(-30 * (0.15 - SEE))) on SSE   ###############################################################

Exp16.GEMSTAT.output.sigmoid_30 <- GEMSTAT_output_Reader_ensemble(.exp.nu = 16, .number.Of.Clusters = 1, copyOutput = F)
Exp16.GEMSTAT.param.perf.sigmoid_30 <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_16/Cluster_1/Outputs", removeDups =F)

aad <- sort(Exp16.GEMSTAT.param.perf.sigmoid_30$performance, decreasing = F, index.return = T)$ix 
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_30, ClusterNumber = 1, paramnumber = 103)
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_30, ClusterNumber = 1, OneGeneManyModel = T, Rank = 1)
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_30, ClusterNumber = 1, OneGeneManyModel = T, GeneNumber = 86)

par(mforw = c(1,1), mar = c(4,4,4,4))
hist(Exp16.GEMSTAT.param.perf.sigmoid_30$performance, breaks = 100)
aa2 <- kmeans(Exp16.GEMSTAT.param.perf.sigmoid_30$parameters, 5)
aa22 <- aa2$cluster
aa <- dist(Exp16.GEMSTAT.param.perf.sigmoid_30$parameters)
aaa <- cmdscale(aa,eig=TRUE, k=2)

# plot solution 
a1 <- aaa$points[,1]
a2 <- aaa$points[,2]
aapch <- rep(0, length(a2))
aapch[aa[1:300]] <- rep(4, 300)
plot(a1, a2, xlab="Coordinate1", ylab="Coordinate2", 
     main="Metric MDS", col = aa22 + 1, pch = aapch)
points(a1[aa[1:300]], a2[aa[1:300]], cex = 3)

######################   EXP 16_using_sigmoid 1/(1 + exp(-3 * (0.15 - SEE))) on SSE   ###############################################################
Exp16.GEMSTAT.output.sigmoid_3 <- GEMSTAT_output_Reader_ensemble(.exp.nu = 16, .number.Of.Clusters = 1, copyOutput = F)
Exp16.GEMSTAT.param.perf.sigmoid_3 <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_16/Cluster_1/Outputs", removeDups =F)

aad <- sort(Exp16.GEMSTAT.param.perf.sigmoid_3$performance, decreasing = F, index.return = T)$ix 
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_3, ClusterNumber = 1, paramnumber = 65)
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_3, ClusterNumber = 1, OneGeneManyModel = T, Rank = 1)
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_3, ClusterNumber = 1, OneGeneManyModel = T, GeneNumber = 61)
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_3, ClusterNumber = 1, OneGeneManyModel = T, plotAmong = c(23,61,65,105,107), varAmong = c(23,61,65,105,107))

Exp16.GEMSTAT.param.perf.sigmoid_3$parameters[c(23,61,65,105,107),]

par(mforw = c(1,1), mar = c(4,4,4,4))
hist(Exp16.GEMSTAT.param.perf.sigmoid_3$performance, breaks = 100)
aa2 <- kmeans(Exp16.GEMSTAT.param.perf.sigmoid_30$parameters, 5)
aa22 <- aa2$cluster
aa <- dist(Exp16.GEMSTAT.param.perf.sigmoid_30$parameters)
aaa <- cmdscale(aa,eig=TRUE, k=2)

# plot solution 
a1 <- aaa$points[,1]
a2 <- aaa$points[,2]
aapch <- rep(0, length(a2))
aapch[aa[1:300]] <- rep(4, 300)
plot(a1, a2, xlab="Coordinate1", ylab="Coordinate2", 
     main="Metric MDS", col = aa22 + 1, pch = aapch)
points(a1[], a2[], cex = 3)
points(a1[65], a2[65], cex = 5)

######################   EXP 16_using_sigmoid 1/(1 + exp(-1 * (0.15 - SEE))) on SSE   ###############################################################
Exp16.GEMSTAT.output.sigmoid_1 <- GEMSTAT_output_Reader_ensemble(.exp.nu = 16, .number.Of.Clusters = 1, copyOutput = F)
Exp16.GEMSTAT.param.perf.sigmoid_1 <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_16/Cluster_1/Outputs", removeDups =F)

aad <- sort(Exp16.GEMSTAT.param.perf.sigmoid_1$performance, decreasing = T, index.return = T)$ix 
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_1, ClusterNumber = 1, paramnumber = 39)
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_1, ClusterNumber = 1, OneGeneManyModel = T, Rank = 1)
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_1, ClusterNumber = 1, OneGeneManyModel = T, GeneNumber = 61)
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_1, ClusterNumber = 1, OneGeneManyModel = T, plotAmong = c(150, 38, 45, 93, 39), varAmong = c(150, 38, 45, 93, 39), Rank = 3)

Exp16.GEMSTAT.param.perf.sigmoid_1$parameters[c(150, 38, 45, 93, 39),]

table(Exp16.GEMSTAT.param.perf.sigmoid_1$parameters[,5])


nrow(Exp16.GEMSTAT.param.perf.sigmoid_1$parameters)

names(Exp16.GEMSTAT.output.sigmoid_1$Allresults[[1]])

rownames(Exp16.GEMSTAT.param.perf.sigmoid_1$parameters)

aa1 <- strsplit(rownames(Exp16.GEMSTAT.param.perf.sigmoid_1$parameters), split = "\\.")
aa2 <- strsplit(names(Exp16.GEMSTAT.output.sigmoid_1$Allresults[[1]]), split = "\\.")

a1 <- unlist(lapply(aa1, '[[', 2))
a2 <- unlist(lapply(aa2, '[[', 2))

######################   EXP 16_using_sigmoid 1/(1 + exp(-1 * (0.15 - SEE))) on SSE   and using exp(- 10 * SSE(enhancer(i))) / sum(exp(- 10 * SSE(enhancer(j)))) as the weight of each enhancer###############################################################

Exp16.GEMSTAT.output.sigmoid_1_softmax_10 <- GEMSTAT_output_Reader_ensemble(.exp.nu = 16, .number.Of.Clusters = 1, copyOutput = F)
Exp16.GEMSTAT.param.perf.sigmoid_1_softmax_10 <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_16/Cluster_1/Outputs", removeDups =F)
Exp16.GEMSTAT.enhWeight.sigmoid_1_softmax_10 <- EnhWeightReader("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_16/Cluster_1/Outputs")

aad <- sort(Exp16.GEMSTAT.param.perf.sigmoid_1_softmax_10$performance, decreasing = T, index.return = T)$ix 
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_1_softmax_10, ClusterNumber = 1, paramnumber = 100)

for(i in 100:length(Exp16.GEMSTAT.enhWeight.sigmoid_1_softmax_10)){
  print(paste("plotting model ", i))
  plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_1_softmax_10, ClusterNumber = 1, paramnumber = i)
  if(! readline("Enter n to go to next model: ") == "n"){
    break
  }
}

plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_1_softmax_10, ClusterNumber = 1, OneGeneManyModel = T, Rank = 1)
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_1_softmax_10, ClusterNumber = 1, OneGeneManyModel = T, GeneNumber = 61)
plot_ensemble_output(Exp16.GEMSTAT.output.sigmoid_1_softmax_10, ClusterNumber = 1, OneGeneManyModel = T, plotAmong = c(150, 38, 45, 93, 39), varAmong = c(150, 38, 45, 93, 39), Rank = 3)

######################   EXP 17_using_sigmoid 1/(1 + exp(-3 * (0.15 - SEE))) on SSE   and using exp(- 5 * SSE(enhancer(i))) / sum(exp(- 5 * SSE(enhancer(j)))) as the weight of each enhancer###############################################################
#adding basal expression to the parameter sampling space
#change the range of parameters in GEMSTAT:
# double ExprPar::min_weight = 0.01;
# double ExprPar::max_weight = 1000;
# double ExprPar::min_interaction = 0.01;	
# double ExprPar::max_interaction = 1000;
# double ExprPar::min_effect_Thermo = 0.01;	
# double ExprPar::max_effect_Thermo = 1000;
# double ExprPar::min_repression = 1.0E-5;
# double ExprPar::max_repression = 1; 
# double ExprPar::min_basal_Thermo = 1.0E-4;	
# double ExprPar::max_basal_Thermo = 1;
#creating range matrix
aa <- matrix(0L, nrow = 7, ncol = 3)
aa[1,] <- c(-2, 1, 1) #ESR1_Binding
aa[2,] <- c(0,  3, 1) #ESR1_txp
aa[3,] <- c(-2, 1, 1) #RARA_binding
aa[4,] <- c(-2, 3, 1) #RARA_txp
aa[5,] <- c(-4, 0, 1) #basal expression
aa[6,] <- c(1,  3, 1) #ESR1_ESR1
aa[7,] <- c(-2, 3, 1) #RARA_ESR1

# creating input for GEMSTAT
GEMSTAT_input_constructor_Ensemble(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2))),
                                   gene_expression_Mat = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2, 
                                   .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                   .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                   .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                   .promoter.sequence=Promoter.Sequence.Char,
                                   .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                   .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                                   .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                   experiment.nu = 17,
                                   .motifs=TF.motifs.TimeSeries.count,
                                   .TF_expression=TF.Expression.TimeSeries,
                                   .max.enh.length = numeric(0),
                                   ..o="DIRECT", coopertingTFIndex=rbind(c(1,1),c(1,2))
                                   , ..oo="SSE", ..ct=12, ..sigma=character(0),
                                   ..na=10
                                   ,.par.RangeMat = aa, .nu.sample=200)


# reading prediction, parameters, and enhancer weights
Exp17.GEMSTAT.output.sigmoid_3_softmax_5 <- GEMSTAT_output_Reader_ensemble(.exp.nu = 17, .number.Of.Clusters = 1, copyOutput = F)
Exp17.GEMSTAT.param.perf.sigmoid_3_softmax_5 <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_17/Cluster_1/Outputs", removeDups =F)
Exp17.GEMSTAT.enhWeight.sigmoid_3_softmax_5 <- EnhWeightReader("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_17/Cluster_1/Outputs")

aad <- sort(Exp17.GEMSTAT.output.sigmoid_3_softmax_5$averageSSE, decreasing = F, index.return = T)$ix 

plot_ensemble_output(Exp17.GEMSTAT.output.sigmoid_3_softmax_5, ClusterNumber = 1, paramnumber = aad[201])
plot_ensemble_output(Exp17.GEMSTAT.output.sigmoid_3_softmax_5, ClusterNumber = 1, OneGeneManyModel = T, Rank = 1)
plot_ensemble_output(Exp17.GEMSTAT.output.sigmoid_3_softmax_5, ClusterNumber = 1, OneGeneManyModel = T, GeneNumber = 61)
plot_ensemble_output(Exp17.GEMSTAT.output.sigmoid_3_softmax_5, ClusterNumber = 1, OneGeneManyModel = T, plotAmong = c(150, 38, 45, 93, 39), varAmong = c(150, 38, 45, 93, 39), Rank = 3)

PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp17.GEMSTAT.output.sigmoid_3_softmax_5,
                    modelIndex = aad[1:19],
                    .Colv = T,
                    .Rowv = T,
                    .dendrogram = "both")
ParameterHeatMap(parameterMat = Exp17.GEMSTAT.param.perf.sigmoid_3_softmax_5$parameters[aad[1:19], ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row")
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp17.GEMSTAT.enhWeight.sigmoid_3_softmax_5,
                     .col_vec=col_vector,
                     .modelIndex=aad[1:19],
                     .removeNonFunctionalRegs=T,
                     .nonFuncThreshold=0.3,
                     .exportplot=T,
                     .dendrogram = "both",
                     .Colv=T,
                     .Rowv = T,
                     .filename = "enhancerHeatmap.png")
Exp17.GEMSTAT.output.sigmoid_3_softmax_5_RegElemWeightCoordinate <- RegElemWeightCoordinate(RegulatoryCoordinateExtractorOutput = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2.Reg.Coor,
                                                                                             EnhWeightReaderOutput = Exp17.GEMSTAT.enhWeight.sigmoid_3_softmax_5, plotheatmap = F)
aaPR2 <- PredictedvsRealCorRMSD(GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.output.sigmoid_3_softmax_5,
                               corrMethod = "pearson" )

aaa2 <-      PlotPercentGenes(..GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.output.sigmoid_3_softmax_5,
                             .LogFileReader.output=Exp17.GEMSTAT.param.perf.sigmoid_3_softmax_5,
                             .EnhWeightReader.output=Exp17.GEMSTAT.enhWeight.sigmoid_3_softmax_5,
                             RMSDrange = c(0,0.15,0.2,1),
                             CorRange = c(1, 0.8, 0.6, -1),
                             .corrMethod = "pearson",
                             Num.Cluster = 5,
                             sortBasedOn = 1,
                             .CorRMSDList = aaPR2,
                             exportplot = T)
#####################   EXP 17_not using sigmoid on SSE   and using exp(- 5 * SSE(enhancer(i))) / sum(exp(- 5 * SSE(enhancer(j)))) as the weight of each enhancer###############################################################
Exp17.GEMSTAT.output.sigmoid_no_softmax_5 <- GEMSTAT_output_Reader_ensemble(.exp.nu = 17, .number.Of.Clusters = 1, copyOutput = F)
Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5 <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_17/Cluster_1/Outputs", removeDups =F)
Exp17.GEMSTAT.enhWeight.sigmoid_no_softmax_5 <- EnhWeightReader("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_17/Cluster_1/Outputs")


aad <- sort(Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5$performance, decreasing = F, index.return = T)$ix 
aad2 <- sort(Exp17.GEMSTAT.output.sigmoid_no_softmax_5$averageSSE, decreasing = F, index.return = T)$ix 

par(mfrow = c(2,1), mar = c(4,4,4,4))
plot(Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5$performance[aad])
abline(h = seq(0.175,0.195, 0.001), col = 2)
abline(v = seq(1,18001, 1000), col = 2)
plot(Exp17.GEMSTAT.output.sigmoid_no_softmax_5$averageSSE[aad2])
abline(h = seq(0.160,0.195, 0.001), col = 2)
abline(v = seq(1,18001, 1000), col = 2)

plot_ensemble_output(Exp17.GEMSTAT.output.sigmoid_no_softmax_5, ClusterNumber = 1, paramnumber = aad2[18200])
plot_ensemble_output(Exp17.GEMSTAT.output.sigmoid_no_softmax_5, ClusterNumber = 1, OneGeneManyModel = T,Rank = 2,groupvsgroup = T,group1 = aad2[1:32],group2 = aad2[15001:15032])

aaa <- RegElemWeightCoordinate(RegulatoryCoordinateExtractorOutput = a, EnhWeightReaderOutput = aa,modelIndex = c(1:19))
ParameterHeatMap(parameterMat = Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5$parameters, exportplot = T, .cellnote = "")
aa <- PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp17.GEMSTAT.output.sigmoid_no_softmax_5,
                    modelIndex = c(aad2[1:100], aad2[16000:16100]), .Colv = F, .Rowv = T, .dendrogram = "row")

GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2.Reg.Coor <- RegulatoryCoordinateExtractor(gene.Names = rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2))
aaa <- RegElemWeightCoordinate(RegulatoryCoordinateExtractorOutput = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2.Reg.Coor,
                               EnhWeightReaderOutput = Exp17.GEMSTAT.enhWeight.sigmoid_no_softmax_5,
                               modelIndex = aad[1:40], removeNonFunctionalRegs = T, nonFuncThreshold = 0.5)

######################Filter the similar models
aaa <- logSimilarityCheck(Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5, logscale = T, thresh = 0.1)
Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5_Filtered <- aaa
Exp17.GEMSTAT.All.Filtered <- FilterForModels(ModelIndex = setdiff(c(1:length(Exp17.GEMSTAT.enhWeight.sigmoid_no_softmax_5)), Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5_Filtered$removed),
                                              GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.output.sigmoid_no_softmax_5,
                                              EnhWeightReader.output = Exp17.GEMSTAT.enhWeight.sigmoid_no_softmax_5,
                                              basedOnSimilarityCheck = F,
                                              .LogFileReader.Output = Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5)
remove(Exp17.GEMSTAT.param.perf.sigmoid_no_softmax_5_Filtered)

aad2 <- sort(Exp17.GEMSTAT.All.Filtered$output$averageSSE, decreasing = F, index.return = T)$ix 
PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp17.GEMSTAT.All.Filtered$output,
                    modelIndex = aad2,
                    .Colv = T,
                    .Rowv = F,
                    .dendrogram = "column")
ParameterHeatMap(parameterMat = Exp17.GEMSTAT.All.Filtered$param$parameters[aad2, ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row")
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp17.GEMSTAT.All.Filtered$enh,
                     .col_vec=col_vector,
                     .modelIndex=aad2,
                     .removeNonFunctionalRegs=F,
                     .nonFuncThreshold=0.3,
                     .exportplot=T,
                     .dendrogram = "both",
                     .Colv=T,
                     .Rowv = T,
                     .filename = "enhancerHeatmap.png")
Exp17.GEMSTAT.All.Filtered_RegElemWeightCoordinate <- RegElemWeightCoordinate(RegulatoryCoordinateExtractorOutput = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2.Reg.Coor,
                                                                              EnhWeightReaderOutput = Exp17.GEMSTAT.All.Filtered$enh, plotheatmap = F)

#MDS of parameters
par(mforw = c(1,1), mar = c(4,4,4,4))
hist(Exp17.GEMSTAT.All.Filtered$param$performance, breaks = 100, main = "Avg Performance (RMSD)")
aa2 <- kmeans(Exp17.GEMSTAT.All.Filtered$param$parameters, 4)
aa22 <- aa2$cluster
aa <- dist(Exp17.GEMSTAT.All.Filtered$param$parameters)
aaa <- cmdscale(aa,eig=TRUE, k=2)

# plot solution 
a1 <- aaa$points[,1]
a2 <- aaa$points[,2]
plot(a1, a2, xlab="Coordinate1", ylab="Coordinate2", 
     main="Metric MDS", col = aa22 + 1)

plot_ensemble_output(Exp17.GEMSTAT.output.sigmoid_no_softmax_5, ClusterNumber = 1, paramnumber = 13)

###Percent gene explained
aaPR <- PredictedvsRealCorRMSD(GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.All.Filtered$output,
                               corrMethod = "pearson" )

aaa <-      PlotPercentGenes(..GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.All.Filtered$output,
                             .LogFileReader.output=Exp17.GEMSTAT.All.Filtered$param,
                             .EnhWeightReader.output=Exp17.GEMSTAT.All.Filtered$enh,
                             RMSDrange = c(0,0.15,0.2,1),
                             CorRange = c(1, 0.8, 0.6, -1),
                             .corrMethod = "pearson",
                             Num.Cluster = 5,
                             sortBasedOn = 1,
                             .CorRMSDList = aaPR,
                             exportplot = T)
Exp17.GEMSTAT.All.Filtered_PercentGenes <- aaa
# since the first cluster of Exp17.GEMSTAT.All.Filtered_PercentGenes seems to be good, I want to
# extract its active enhancers and see if they are transcribed in GRO-Seq.R script
aa <- which(Exp17.GEMSTAT.All.Filtered_PercentGenes$Clustering$enhBased$cluster %in% 1)
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp17.GEMSTAT.All.Filtered$enh,
                     .col_vec=col_vector,
                     .modelIndex=aa,
                     .removeNonFunctionalRegs=F,
                     .nonFuncThreshold=0.3,
                     .exportplot=T,
                     .dendrogram = "column",
                     .Colv=T,
                     .Rowv = F,
                     .filename = "enhancerHeatmap.png")
ParameterHeatMap(parameterMat = Exp17.GEMSTAT.All.Filtered$param$parameters[aa, ],
                 exportplot = T,
                 .Rowv = F,
                 .dendrogram = "none")

PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp17.GEMSTAT.All.Filtered$output,
                    modelIndex = aa,
                    .Colv = T,
                    .Rowv = T,
                    .dendrogram = "both")
#
aa <- sort(Exp17.GEMSTAT.All.Filtered_PercentGenes$PercentExplained$Cor[1,], decreasing = T, index.return = T)$ix[1:12]

plot_ensemble_output(Exp17.GEMSTAT.All.Filtered$output, ClusterNumber = 1, paramnumber = 977)

aaa <- Exp17.GEMSTAT.All.Filtered_RegElemWeightCoordinate[, c(c(1:5), (aa + 5))]
aaaa <- sort(rowSums(Exp17.GEMSTAT.All.Filtered_RegElemWeightCoordinate[,  aa + 5]), decreasing = T, index.return= T)
par(mfrow =c(1,1), mar = c(4,4,4,4))
par(xpd=FALSE)
plot(aaaa$x)
abline(v=c(60, 80, 100,120,140,160), col = 2)
aaaaa <- aaa[aaaa$ix[1:200], c(1:5)]
aaaaa <- aaaaa[which(aaaaa[,4]!=1001), ]
aaaaa.GR <- makeGRangesFromDataFrame(aaaaa)
aaaaaa <- overlapExtractor(aaaaa.GR, GSE67295.GROseq.GRanges)
unique(unlist(aaaaaa))
aaaaaaa <- GSE67295.GROseq[unique(unlist(aaaaaa)), ]
boxplot.matrix(as.matrix(aaaaaaa[,9:16]), las = 2, outline = F, ylim = c(0,5000))
boxplot.matrix(as.matrix(GSE67295.GROseq[,9:16]), las = 2, outline = F)

aaaaaa2 <- overlapExtractor(aaaaa.GR, GSE27463.GROseq.GRanges)
unique(unlist(aaaaaa2))
aaaaaaa2 <- GSE27463.GROseq[unique(unlist(aaaaaa2)),]
boxplot.matrix(as.matrix(aaaaaaa2[,9:16]), las = 2)
##########################################################################################
#sort models by performance:
aa <- sort(Exp17.GEMSTAT.All.Filtered$output$averageSSE, decreasing = F, index.return = T)$ix
aa <- aa[1:100]
aaa <- (Exp17.GEMSTAT.All.Filtered_RegElemWeightCoordinate[,  aa + 5])
aaaa <- numeric(nrow(aaa))
for(i in 1:nrow(aaa)){
  aaaa[i] <- median(as.numeric(aaa[i, ]))
}
#take 27 positive and 27 negative examples
aaPos <- Exp17.GEMSTAT.All.Filtered_RegElemWeightCoordinate[which(aaaa > 0.5), 1:5]
aaNeg <- Exp17.GEMSTAT.All.Filtered_RegElemWeightCoordinate[which(aaaa < 0.01)[1:27], 1:5]

#look at their transcription
aaPos.GR <- makeGRangesFromDataFrame(aaPos)
aaNeg.GR <- makeGRangesFromDataFrame(aaNeg)
aaPos.Overlap1 <- overlapExtractor(aaPos.GR, GSE67295.GROseq.GRanges)
aaPos.Overlap2 <- overlapExtractor(aaPos.GR, GSE27463.GROseq.GRanges)
aaNeg.Overlap1 <- overlapExtractor(aaNeg.GR, GSE67295.GROseq.GRanges)
aaNeg.Overlap2 <- overlapExtractor(aaNeg.GR, GSE27463.GROseq.GRanges)

par(mfrow = c(1,2), mar = c(4,4,1,1))
# aaaaaaa <- as.matrix(GSE67295.GROseq[unique(unlist(aaPos.Overlap1)),c(9:16)])
# boxplot.matrix(aaaaaaa, las = 2, outline = F , xaxt = "n", ylim = c(1,2500))
# aaaaaaa <- as.matrix(GSE67295.GROseq[unique(unlist(aaNeg.Overlap1)),c(9:16)], xaxt = "n")
# boxplot.matrix(aaaaaaa, las = 2, outline = F, xaxt = "n", ylim = c(1,2500))
aaaaaaa <- as.matrix(GSE27463.GROseq[unique(unlist(aaPos.Overlap2)),c(7:10)])
boxplot.matrix(aaaaaaa, las = 2, outline = F, xaxt = "n", ylim = c(1,200))
aaaaaaa <- as.matrix(GSE27463.GROseq[unique(unlist(aaNeg.Overlap2)),c(7:10)])
boxplot.matrix(aaaaaaa, las = 2, outline = F, xaxt = "n", ylim = c(1,200))



unique(unlist(aaNeg.Overlap1))

aaaaaaa <- as.matrix(GSE67295.GROseq[unique(unlist(aaPos.Overlap1)),c(9:16)])
boxplot.matrix(as.matrix(aaaaaaa[,9:16]), las = 2, outline = F, ylim = c(0,5000))
boxplot.matrix(as.matrix(GSE67295.GROseq[,9:16]), las = 2, outline = F)

#########
#Filter the models which have ER as a repressor or have less than one coop between ER-ER

aa <- union(which(Exp17.GEMSTAT.All.Filtered$param$parameters[,2] < 1), which(Exp17.GEMSTAT.All.Filtered$param$parameters[,6] < 1))
Exp17.GEMSTAT.All.Filtered.NoERrepression <- FilterForModels(ModelIndex = setdiff(c(1:length(Exp17.GEMSTAT.All.Filtered$enh)), aa),
                                              GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.All.Filtered$output,
                                              EnhWeightReader.output = Exp17.GEMSTAT.All.Filtered$enh,
                                              basedOnSimilarityCheck = F,
                                              .LogFileReader.Output = Exp17.GEMSTAT.All.Filtered$param)

aad2 <- sort(Exp17.GEMSTAT.All.Filtered.NoERrepression$output$averageSSE, decreasing = F, index.return = T)$ix 
PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp17.GEMSTAT.All.Filtered.NoERrepression$output,
                    modelIndex = aad2,
                    .Colv = F,
                    .Rowv = F,
                    .dendrogram = "none")
ParameterHeatMap(parameterMat = Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters[aad2, ],
                 exportplot = T,
                 .Rowv = F,
                 .dendrogram = "none")
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp17.GEMSTAT.All.Filtered.NoERrepression$enh,
                     .col_vec=col_vector,
                     .modelIndex=aad2,
                     .removeNonFunctionalRegs=T,
                     .nonFuncThreshold=0.5,
                     .exportplot=T,
                     .dendrogram = "row",
                     .Colv=F,
                     .Rowv = T,
                     .filename = "enhancerHeatmap.png")
Exp17.GEMSTAT.All.Filtered.NoERrepression_RegElemWeightCoordinate <- RegElemWeightCoordinate(RegulatoryCoordinateExtractorOutput = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2.Reg.Coor,
                                                                              EnhWeightReaderOutput = Exp17.GEMSTAT.All.Filtered.NoERrepression$enh, plotheatmap = F)

#MDS of parameters
par(mforw = c(1,1), mar = c(4,4,4,4))
hist(Exp17.GEMSTAT.All.Filtered.NoERrepression$param$performance, breaks = 100, main = "Avg Performance (RMSD)")
aa2 <- kmeans(Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters, 4)
aa22 <- aa2$cluster
aa <- dist(Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters)
aaa <- cmdscale(aa,eig=TRUE, k=2)

# plot solution 
a1 <- aaa$points[,1]
a2 <- aaa$points[,2]
plot(a1, a2, xlab="Coordinate1", ylab="Coordinate2", 
     main="Metric MDS", col = aa22 + 1)

plot_ensemble_output(Exp17.GEMSTAT.All.Filtered.NoERrepression$output, ClusterNumber = 1, paramnumber = 513)

###Percent gene explained
aaPR <- PredictedvsRealCorRMSD(GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.All.Filtered.NoERrepression$output,
                               corrMethod = "pearson" )

aaa <-      PlotPercentGenes(..GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.All.Filtered.NoERrepression$output,
                             .LogFileReader.output=Exp17.GEMSTAT.All.Filtered.NoERrepression$param,
                             .EnhWeightReader.output=Exp17.GEMSTAT.All.Filtered.NoERrepression$enh,
                             RMSDrange = c(0,0.15,0.2,1),
                             CorRange = c(1, 0.8, 0.6, -1),
                             .corrMethod = "pearson",
                             Num.Cluster = 5,
                             sortBasedOn = 1,
                             .CorRMSDList = aaPR,
                             exportplot = T)
Exp17.GEMSTAT.All.Filtered.NoERrepression_PercentGenes <- aaa

### cluster 3 of enhancer based in Exp17.GEMSTAT.All.Filtered.NoERrepression_PercentGenes seems to be better than others
aa <- which(Exp17.GEMSTAT.All.Filtered.NoERrepression_PercentGenes$Clustering$enhBased$cluster %in% 3)
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp17.GEMSTAT.All.Filtered.NoERrepression$enh,
                     .col_vec=col_vector,
                     .modelIndex=aa,
                     .removeNonFunctionalRegs=T,
                     .nonFuncThreshold=0.3,
                     .exportplot=T,
                     .dendrogram = "none",
                     .Colv=F,
                     .Rowv = T,
                     .filename = "enhancerHeatmap.png")
ParameterHeatMap(parameterMat = Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters[aa, ],
                 exportplot = T,
                 .Rowv = F,
                 .dendrogram = "none")

PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp17.GEMSTAT.All.Filtered.NoERrepression$output,
                    modelIndex = aa,
                    .Colv = F,
                    .Rowv = T,
                    .dendrogram = "none")

#sort based on the last parameter : ER-RAR and look at the models again
aasp7 <- sort(Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters[aa,7], decreasing = T, index.return = T)$ix

RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp17.GEMSTAT.All.Filtered.NoERrepression$enh,
                     .col_vec=col_vector,
                     .modelIndex=aa[aasp7],
                     .removeNonFunctionalRegs=T,
                     .nonFuncThreshold=0.3,
                     .exportplot=T,
                     .dendrogram = "none",
                     .Colv=F,
                     .Rowv = F,
                     .filename = "enhancerHeatmap.png")
ParameterHeatMap(parameterMat = Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters[aa[aasp7], ],
                 exportplot = T,
                 .Rowv = F,
                 .dendrogram = "none")

PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp17.GEMSTAT.All.Filtered.NoERrepression$output,
                    modelIndex = aa[aasp7],
                    .Colv = F,
                    .Rowv = F,
                    .dendrogram = "none")

plot_ensemble_output(Exp17.GEMSTAT.All.Filtered.NoERrepression$output, ClusterNumber = 1, OneGeneManyModel = T,varAmong = aa, plotAmong = aa)


CompareModelPlot(GEMSTAT_output_Reader_ensemble.Output = Exp17.GEMSTAT.All.Filtered.NoERrepression$output,
                 LogFileReader.Output = Exp17.GEMSTAT.All.Filtered.NoERrepression$param,
                 EnhWeightReader.Output = Exp17.GEMSTAT.All.Filtered.NoERrepression$enh,
                 modelIndex = c(295,802),
                 ..col_vec = col_vector,
                 .geneNames=rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2),
                 ..exportplot = T,
                 plotDir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Ensemble/EXP17")


#######find groups of models which are only different in ER-RAR coop:
a1 <- which(Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters[,1] < 0.1)
a2 <- which(Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters[,2] > 1)
a3 <- which(Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters[,3] < 0.1)
a4 <- which(Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters[,4] > 1)
a5 <- which(Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters[,6] > 1)
a6 <- which(Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters[,7] > 1)
aaRARup <- Reduce(intersect, list(a1, a2, a3, a4 , a5, a6))

a6 <- which(Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters[,7] < 1)
aaRARdown <- Reduce(intersect, list(a1, a2, a3, a4 , a5, a6))
par(mfrow = c(2,1), mar = c(4,4,4,4))
hist(Exp17.GEMSTAT.All.Filtered.NoERrepression$output$averageSSE[aaRARup], breaks = 50)
hist(Exp17.GEMSTAT.All.Filtered.NoERrepression$output$averageSSE[aaRARdown], breaks = 50)

aa <- which(Exp17.GEMSTAT.All.Filtered.NoERrepression$output$averageSSE[aaRARup] < 0.17)
aaRARup <- aaRARup[aa]
aa <- which(Exp17.GEMSTAT.All.Filtered.NoERrepression$output$averageSSE[aaRARdown] < 0.17)
aaRARdown <- aaRARdown[aa]


RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp17.GEMSTAT.All.Filtered.NoERrepression$enh,
                     .col_vec=col_vector,
                     .modelIndex=c(aaRARup,aaRARdown),
                     .removeNonFunctionalRegs=T,
                     .nonFuncThreshold=0.3,
                     .exportplot=T,
                     .dendrogram = "both",
                     .Colv=T,
                     .Rowv = T,
                     .filename = "enhancerHeatmap.png")
ParameterHeatMap(parameterMat = Exp17.GEMSTAT.All.Filtered.NoERrepression$param$parameters[c(aaRARup,aaRARdown), ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row")

PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp17.GEMSTAT.All.Filtered.NoERrepression$output,
                    modelIndex = c(aaRARup,aaRARdown),
                    .Colv = F,
                    .Rowv = T,
                    .dendrogram = "row")

plot_ensemble_output(Exp17.GEMSTAT.All.Filtered.NoERrepression$output, ClusterNumber = 1, OneGeneManyModel = T,groupvsgroup = T, group1 = aaRARup, group2 = aaRARdown, Rank = 1)


CompareModelPlot(GEMSTAT_output_Reader_ensemble.Output = Exp17.GEMSTAT.All.Filtered.NoERrepression$output,
                 LogFileReader.Output = Exp17.GEMSTAT.All.Filtered.NoERrepression$param,
                 EnhWeightReader.Output = Exp17.GEMSTAT.All.Filtered.NoERrepression$enh,
                 modelIndex = c(3,2),
                 ..col_vec = col_vector,
                 .geneNames=rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2),
                 ..exportplot = T,
                 plotDir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Ensemble/EXP17")

# models 2 and 3 can be compared on genes: 7057, 7494, 8204, 7779, 91107, 333926

################################################################################################################################################################################
################################################################################################################################################################################
#####################   EXP 17_not using sigmoid on SSE   and using -ln(sum(exp(- 5 * SSE(enhancer(j))))) as the obj function###############################################################
Exp17.GEMSTAT.output.sigmoid_no_objss_5 <- GEMSTAT_output_Reader_ensemble(.exp.nu = 17, .number.Of.Clusters = 1, copyOutput = F)
Exp17.GEMSTAT.param.perf.sigmoid_no_objss_5 <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_17/Cluster_1/Outputs", removeDups =F)
Exp17.GEMSTAT.enhWeight.sigmoid_no_objss_5 <- EnhWeightReader("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_17/Cluster_1/Outputs")


aad <- sort(Exp17.GEMSTAT.output.sigmoid_no_objss_5$averageSSE, decreasing = F, index.return = T)$ix 

par(mfrow = c(1,1), mar = c(4,4,4,4))
plot(Exp17.GEMSTAT.output.sigmoid_no_objss_5$averageSSE[aad], main = "obj_ss_5_performance", ylab = "SSE", xlab="models")
abline(h = seq(0.150,0.195, 0.001), col = 2)
abline(v = seq(1,18001, 1000), col = 2)

plot_ensemble_output(Exp17.GEMSTAT.output.sigmoid_no_objss_5, ClusterNumber = 1, paramnumber = aad[1])

######################Filter the similar models
aaa <- logSimilarityCheck(Exp17.GEMSTAT.param.perf.sigmoid_no_objss_5, logscale = T, thresh = 0.1)
Exp17.GEMSTAT.param.perf.sigmoid_no_objss_5_Filtered <- aaa
Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5 <- FilterForModels(ModelIndex = setdiff(c(1:length(Exp17.GEMSTAT.output.sigmoid_no_objss_5$averageSSE)), Exp17.GEMSTAT.param.perf.sigmoid_no_objss_5_Filtered$removed),
                                              GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.output.sigmoid_no_objss_5,
                                              EnhWeightReader.output = Exp17.GEMSTAT.enhWeight.sigmoid_no_objss_5,
                                              basedOnSimilarityCheck = F,
                                              .LogFileReader.Output = Exp17.GEMSTAT.param.perf.sigmoid_no_objss_5)
remove(Exp17.GEMSTAT.param.perf.sigmoid_no_objss_5_Filtered)

aad <- sort(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output$averageSSE, decreasing = F, index.return = T)$ix 

par(mfrow = c(1,1), mar = c(4,4,4,4))
plot(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output$averageSSE[aad], main = "obj_ss_5_performance_Filtered", ylab = "SSE", xlab="models")
abline(h = seq(0.150,0.195, 0.001), col = 2)
#########

PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output,
                    modelIndex = aad,
                    .Colv = F,
                    .Rowv = F,
                    .dendrogram = "none")
ParameterHeatMap(parameterMat = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[aad, ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row")
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$enh,
                     .removeNonFunctionalRegs = T, .nonFuncThreshold = 0.1,
                     .exportplot = T,
                     .modelIndex = aad,
                     .dendrogram = "none",
                     .Colv = F, .Rowv = F,
                     setmaxto1RestZero = T)

Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5_RegElemWeightCoordinate <- RegElemWeightCoordinate(RegulatoryCoordinateExtractorOutput = GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2.Reg.Coor,
                                                                              EnhWeightReaderOutput = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$enh, plotheatmap = F, setmaxto1RestZero = T)

###Percent gene explained
aaPR <- PredictedvsRealCorRMSD(GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output,
                               corrMethod = "pearson" )

aaa <-      PlotPercentGenes(..GEMSTAT_output_Reader_ensemble.output = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output,
                             .LogFileReader.output=Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param,
                             .EnhWeightReader.output=Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$enh,
                             RMSDrange = c(0,0.15,0.2,1),
                             CorRange = c(1, 0.8, 0.6, -1),
                             .corrMethod = "pearson",
                             Num.Cluster = 4,
                             sortBasedOn = 1,
                             .CorRMSDList = aaPR,
                             exportplot = T)

###########################
#######find groups of models which are only different in ER-RAR coop:
a1 <- which(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[,1] < 0.1)
a2 <- which(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[,2] > 1)
a3 <- which(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[,3] < 0.1)
a4 <- which(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[,4] > 1)
a5 <- which(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[,6] > 1)
a6 <- which(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[,7] > 10)
aaRARup <- Reduce(intersect, list(a1, a2, a3, a4 , a5, a6))

a6 <- which(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[,7] < 0.1)
aaRARdown <- Reduce(intersect, list(a1, a2, a3, a4 , a5, a6))
par(mfrow = c(2,1), mar = c(4,4,4,4))
hist(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output$averageSSE[aaRARup], breaks = 50)
hist(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output$averageSSE[aaRARdown], breaks = 50)

aa <- which(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output$averageSSE[aaRARup] < 0.17)
aaRARup <- aaRARup[aa]
aa <- which(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output$averageSSE[aaRARdown] < 0.17)
aaRARdown <- aaRARdown[aa]

aaacolside <- c(rep("blue", length(aaRARup)), rep("pink", length(aaRARdown)))
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$enh,
                     .col_vec=col_vector,
                     .modelIndex=c(aaRARup,aaRARdown),
                     .removeNonFunctionalRegs=T,
                     .nonFuncThreshold=0.3,
                     .exportplot=T,
                     .dendrogram = "row",
                     .Colv=F,
                     .Rowv = T,
                     setmaxto1RestZero = T,
                     .filename = "enhancerHeatmap.png",
                     .ColSideColors = aaacolside)
ParameterHeatMap(parameterMat = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[c(aaRARup,aaRARdown), ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row",
                 .RowSideColors = aaacolside)

PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output,
                    modelIndex = c(aaRARup,aaRARdown),
                    .Colv = F,
                    .Rowv = F,
                    .dendrogram = "none",
                    .ColSideColors = aaacolside)

plot_ensemble_output(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output, ClusterNumber = 1, OneGeneManyModel = T,groupvsgroup = T, group1 = aaRARup, group2 = aaRARdown, Rank = 1)


#genes 56 and 45 can be good examples

names(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$enh)[aaRARdown]
names(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$enh)[aaRARup]

CompareModelPlot(GEMSTAT_output_Reader_ensemble.Output = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output,
                 LogFileReader.Output = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param,
                 EnhWeightReader.Output = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$enh,
                 modelIndex = c(14,117),
                 ..col_vec = col_vector,
                 .geneNames=rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2),
                 ..exportplot = T,
                 plotDir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Ensemble/EXP17",
                 .setmaxto1RestZero = T)
#See if for the genes that are explained well: the chosen enhancer changes from model to model : models to compare: 14, 117
# take the parameters of these two models and run GEMSTAT in multiple different purterbation settings: e.g : setting the binding of each TF to zero, setting txp to 1, coop to 1 and ...

PurturbationInputCreator(modelParameters = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[14,],
                         TFnames=c("ESR1", "RARA"),
                         .coopertingTFIndex= rbind(c(1,1),c(1,2)),
                         modelIndex=14,
                         experiment.nu = 17,
                         .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                         .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                         .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/")
PurturbationInputCreator(modelParameters = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[117,],
                         TFnames=c("ESR1", "RARA"),
                         .coopertingTFIndex= rbind(c(1,1),c(1,2)),
                         modelIndex=117,
                         experiment.nu = 17,
                         .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                         .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_TrySoftMax/seq2expr",
                         .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/")

Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5_output_Purt_14_117 <- Purturbation_output_Reader (.exp.nu = 17,
                                   modelIndex = 14,
                                   .root_dir="tabebor2@veda.cs.illinois.edu:/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble",
                                   .dest_dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble"
                                   ,copyOutput = T)
Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5_Parameters_Purt_14_117 <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_17/Cluster_1/Outputs/PurturbOutputs", removeDups =F)
Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5_EnhWeight_Purt_14_117 <- EnhWeightReader("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_17/Cluster_1/Outputs/PurturbOutputs")

plot_ensemble_output(Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5_output_Purt_14_117, ClusterNumber = 1, paramnumber = 10)

aaacolside <- c(rep("blue", 5), rep("pink", 5))
PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5_output_Purt_14_117,
                    .Colv = F,
                    .Rowv = T,
                    .dendrogram = "row",
                    .ColSideColors = aaacolside)
RegElemWeightHeatMap(.EnhWeightReaderOutput =Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5_EnhWeight_Purt_14_117,
                     .col_vec=col_vector,
                     #.modelIndex=c(aaRARup,aaRARdown),
                     .removeNonFunctionalRegs=F,
                     .nonFuncThreshold=0.3,
                     .exportplot=T,
                     .dendrogram = "row",
                     .Colv=F,
                     .Rowv = T,
                     setmaxto1RestZero = T,
                     .filename = "enhancerHeatmap.png",
                     .ColSideColors = aaacolside)
ParameterHeatMap(parameterMat = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5_Parameters_Purt_14_117$parameters,
                 exportplot = T,
                 .Rowv = F,
                 .dendrogram = "none",
                 .RowSideColors = aaacolside)

CompareModelPlot(GEMSTAT_output_Reader_ensemble.Output = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5_output_Purt_14_117,
                 LogFileReader.Output = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5_Parameters_Purt_14_117,
                 EnhWeightReader.Output = Exp17.GEMSTAT.All.Filtered.sigmoid_no_objss_5_EnhWeight_Purt_14_117,
                 modelIndex = c(6,10),
                 ..col_vec = col_vector,
                 .geneNames=rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2),
                 ..exportplot = T,
                 plotDir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Ensemble/EXP17",
                 .setmaxto1RestZero = T)
########################################################################################################################
########################################################################################################################
#Seems like the cooperation parameter is not really being used
#Run the ensemble with new motifs: dimer motifs for both

##################################################################################
#create motif count matrix for ESR1 and RARA : use dimer sites for both of them
TF.motifs.TimeSeries.Dimer <- list()
aa <- list.files(path="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Motifs/Chosen_timeSeries_dimer/", full.names=T, recursive=FALSE)
for( i in 1:length(aa)){
  TF.motifs.TimeSeries.Dimer[[i]] <- read.table(aa[i], header=T, sep ="\t")
  TF.motifs.TimeSeries.Dimer[[i]] =  TF.motifs.TimeSeries.Dimer[[i]][,2:5]
}

names(TF.motifs.TimeSeries.Dimer) <- c("ESR1","RARA")
TF.motifs.TimeSeries.Dimer.count <- lapply(TF.motifs.TimeSeries.Dimer, PWMtoCount)
#########
#write the motifs to a file:
MotifWriter(motif.List = TF.motifs.TimeSeries.Dimer.count, output.File.Name = "motifs_timeseries_dimer")
##########
## Left this for later
########################################################################################################################
########################################################################################################################
######################   EXP 18_using -ln(sum(exp(- 5 * SSE(enhancer(j))))) as the obj functionr###############################################################
#using genes found in motif analysis
ER.associated.genes.MotifAnalysis.Timeseries.Exp.Mat
aa1 <- match(ER.associated.genes.ERERgt2.ERRAReq0.GEMSTAT, rownames(GSE78167.RNAseq.Avg.Norm))
aa2 <- match(ER.associated.genes.ERERgt3.ERRARgt3.GEMSTAT, rownames(GSE78167.RNAseq.Avg.Norm))
aa1 <- aa1[!is.na(aa1)]
aa2 <- aa2[!is.na(aa2)]

sum(ER.associated.genes.ERERgt2.ERRAReq0.GEMSTAT %in% rownames(GSE78167.RNAseq.Diff))
sum(ER.associated.genes.ERERgt3.ERRARgt3.GEMSTAT %in% rownames(GSE78167.RNAseq.Diff))

GSE78167.RNAseq.Diff

par(mfrow = c(5, 5), mar= c(0.1, 0.1, 0.1, 0.1))
for(i in 1:length(aa1)){
  plot(GSE78167.RNAseq.Avg.Norm[aa1[i], ], type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main="")
}
par(mfrow = c(5, 5), mar= c(0.1, 0.1, 0.1, 0.1))
for(i in 1:length(aa2)){
  plot(GSE78167.RNAseq.Avg.Norm[aa2[i], ], type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main="")
}
ER.associated.genes.MotifBased <- c(ER.associated.genes.ERERgt2.ERRAReq0.GEMSTAT, ER.associated.genes.ERERgt3.ERRARgt3.GEMSTAT)
#creating the expression matrix:
aa1 <- match(ER.associated.genes.MotifBased, rownames(GSE78167.RNAseq.Avg.Norm))
ER.associated.genes.MotifBased.ExpMat <- GSE78167.RNAseq.Avg.Norm[aa1, ]

#creating the motifs:
TF.motifs.TimeSeries.loose <- lapply(X = TF.motifs.TimeSeries.loose.t, t) 
TF.motifs.TimeSeries.loose.count <- lapply(TF.motifs.TimeSeries.loose, PWMtoCount)

#adding basal expression to the parameter sampling space
#change the range of parameters in GEMSTAT:
# double ExprPar::min_weight = 0.01;
# double ExprPar::max_weight = 1000;
# double ExprPar::min_interaction = 0.01;	
# double ExprPar::max_interaction = 1000;
# double ExprPar::min_effect_Thermo = 0.01;	
# double ExprPar::max_effect_Thermo = 1000;
# double ExprPar::min_repression = 1.0E-5;
# double ExprPar::max_repression = 1; 
# double ExprPar::min_basal_Thermo = 1.0E-4;	
# double ExprPar::max_basal_Thermo = 1;
#creating range matrix
aa <- matrix(0L, nrow = 7, ncol = 3)
aa[1,] <- c(-2, 1, 1) #ESR1_Binding
aa[2,] <- c(0,  3, 1) #ESR1_txp
aa[3,] <- c(-2, 1, 1) #RARA_binding
aa[4,] <- c(-2, 3, 1) #RARA_txp
aa[5,] <- c(-4, 0, 1) #basal expression
aa[6,] <- c(1,  3, 1) #ESR1_ESR1
aa[7,] <- c(-2, 3, 1) #RARA_ESR1

#ex
# creating input for GEMSTAT
GEMSTAT_input_constructor_Ensemble(clusterIndex = list(c(1:nrow(ER.associated.genes.MotifBased.ExpMat))),
                                   gene_expression_Mat = ER.associated.genes.MotifBased.ExpMat, 
                                   .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                   .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                   .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                   .promoter.sequence=Promoter.Sequence.Char,
                                   .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                   .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_Max/seq2expr",
                                   .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                   experiment.nu = 18,
                                   .par.exp.nu = 17,
                                   .motifs=TF.motifs.TimeSeries.loose.count,
                                   .TF_expression=TF.Expression.TimeSeries,
                                   .max.enh.length = numeric(0),
                                   ..o="DIRECT", coopertingTFIndex=rbind(c(1,1),c(1,2))
                                   , ..oo="SSE", ..ct=11, ..sigma=character(0),
                                   ..na=10
                                   ,.par.RangeMat = aa,
                                   .nu.sample=200,
                                   writeParameters = F)


# reading prediction, parameters, and enhancer weights
Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5 <- list()
Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]] <- GEMSTAT_output_Reader_ensemble(.exp.nu = 18, .number.Of.Clusters = 1, copyOutput = F)
Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]] <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_18/Cluster_1/Outputs", removeDups =F)
Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[3]] <- EnhWeightReader("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_18/Cluster_1/Outputs")
names(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5) <- c("output", "param", "enh")

aad <- sort(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE, decreasing = F, index.return = T)$ix 

plot_ensemble_output(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                     ClusterNumber = 1,
                     paramnumber = aad[1])
plot_ensemble_output(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                     ClusterNumber = 1,
                     paramnumber = 278)
plot_ensemble_output(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                     ClusterNumber = 1,
                     OneGeneManyModel = T, Rank = 1)
plot_ensemble_output(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                     ClusterNumber = 1,
                     OneGeneManyModel = T, GeneNumber = 61)
plot_ensemble_output(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                     ClusterNumber = 1,
                     OneGeneManyModel = T,
                     plotAmong = c(150, 38, 45, 93, 39), varAmong = c(150, 38, 45, 93, 39), Rank = 3)

######################Filter the similar models
aaa <- logSimilarityCheck(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]], logscale = T, thresh = 0.1)
Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5_filterpar <- aaa
Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5 <- FilterForModels(ModelIndex = setdiff(c(1:length(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE)), Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5_filterpar$removed),
                                                                 GEMSTAT_output_Reader_ensemble.output = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                                                                 EnhWeightReader.output = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[3]],
                                                                 basedOnSimilarityCheck = F,
                                                                 .LogFileReader.Output = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]])
remove(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5_filterpar)
aad <- sort(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE, decreasing = F, index.return = T)$ix 

par(mfrow = c(1,1), mar = c(4,4,4,4))
plot(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE[aad], main = "obj_ss_5_performance_Filtered", ylab = "SSE", xlab="models")
abline(h = seq(0.160,0.230, 0.001), col = "grey", lty = 3, lwd = 0.7)
#########
aaarowsidecol <- c(rep("blue", 22), rep("pink", 22))
PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                    modelIndex = aad,
                    .Colv = F,
                    .Rowv = F,
                    .dendrogram = "none",
                    .RowSideColors = aaarowsidecol)
aadd <- which(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE < 0.18)
ParameterHeatMap(parameterMat = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[aadd, ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row")
###################################################
#create the charachter vector of the sequence of all reg elements of the modeled genes:
aaa <- match(rownames(ER.associated.genes.MotifBased.ExpMat), names(Vicinity100kb.Enhancer.plusPromoter.By.gene))
ER.associated.genes.MotifBased.reg.elements.seqANDlength <- WriteFastaOfBag(enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene[aaa],
                                              Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                              promoter.promoter.int=Promoter.gene.PromoterIntList[aaa],
                                              promoter.sequence=Promoter.Sequence.Char,
                                              return.Length = T,
                                              returnListNotwrite = T)
ER.associated.genes.MotifBased.reg.elements.seq <- do.call(c, ER.associated.genes.MotifBased.reg.elements.seq$sequence)

##sequence of regulatory elements of interest:
ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence
ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence
# ER.associated.genes.MotifBased.reg.elements.Interest is an integer vector with one entry associated with each reg element. entry is 1: ER no RAR, entry is 2: ER and RAR, otherwise zero
ER.associated.genes.MotifBased.reg.elements.Interest <- integer(length = length(ER.associated.genes.MotifBased.reg.elements.seq))
#Find which of the reg elements are the ones of interest, assign color blue to ones without RAR and pink to ones with RAR, and white to others
aarowsidecol <- rep("white", length(ER.associated.genes.MotifBased.reg.elements.seq))
aa <- which(toupper(ER.associated.genes.MotifBased.reg.elements.seq) %in% toupper(ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence))
aarowsidecol[aa] <- "pink"
ER.associated.genes.MotifBased.reg.elements.Interest[aa] <- 2
aa <- which(toupper(ER.associated.genes.MotifBased.reg.elements.seq) %in% toupper(ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence))
aarowsidecol[aa] <- "blue"
ER.associated.genes.MotifBased.reg.elements.Interest[aa] <- 1
###################################################################
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[3]],
                     .removeNonFunctionalRegs = T, .nonFuncThreshold = 0.1,
                     .exportplot = T,
                     .modelIndex = aad,
                     .dendrogram = "column",
                     .Colv = T, .Rowv = F,
                     setmaxto1RestZero = T
                     , .RowSideCol = aarowsidecol
                     )
GSE78167.RNAseq.Avg.ChosenGenes.Norm.Reg.MotifBased.Coor <- RegulatoryCoordinateExtractor(gene.Names = rownames(ER.associated.genes.MotifBased.ExpMat))
Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5_RegElemWeightCoordinate <- RegElemWeightCoordinate(RegulatoryCoordinateExtractorOutput = GSE78167.RNAseq.Avg.ChosenGenes.Norm.Reg.MotifBased.Coor,
                                                                                                 EnhWeightReaderOutput = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[3]],
                                                                                                 plotheatmap = F,
                                                                                                 setmaxto1RestZero = T)
####Look which models use the enhancers of interest the most: THink how to plot this??
# for each gene compare the average performance of models which use enhancer of interest vs other enhancers
aa <- integer(length(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[3]]))
for( i in 1: length(aa)){
  aa[i] <- sum(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5_RegElemWeightCoordinate[[5 + i]][which(ER.associated.genes.MotifBased.reg.elements.Interest > 0)])
}
hist(aa ,main = "no. of enhancers of interest used/out of 61")

###Percent gene explained
aaPR <- PredictedvsRealCorRMSD(GEMSTAT_output_Reader_ensemble.output = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                               corrMethod = "pearson" )

aaa <-      PlotPercentGenes(..GEMSTAT_output_Reader_ensemble.output = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                             .LogFileReader.output=Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]],
                             .EnhWeightReader.output=Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[3]],
                             RMSDrange = c(0,0.15,0.2,1),
                             CorRange = c(1, 0.8, 0.6, -1),
                             .corrMethod = "pearson",
                             Num.Cluster = 4,
                             sortBasedOn = 1,
                             .CorRMSDList = aaPR,
                             exportplot = T)

#find models with ER as activator, ER-ER interaction positive
aaERact <- which(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[, 2] > 1)
aaERERpos <- which(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[, 6] > 1)
aaERRARpos <- which(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[, 7] > 1)
aaERRARneg <- which(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[, 7] < 1)

aaERactERERpos <- intersect(aaERact, aaERERpos)
aaERactERERposERRARpos <- intersect(aaERactERERpos, aaERRARpos)
aaERactERERposERRARneg <- intersect(aaERactERERpos, aaERRARneg)
hist(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE[aaERactERERpos])
aadd <- which(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE < 0.18)

ParameterHeatMap(parameterMat = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[aaERactERERpos, ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row")

aaacolside <- c(rep("pink", length(aaERactERERposERRARpos)), rep("blue", length(aaERactERERposERRARneg)))
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5$enh,
                     .col_vec=col_vector,
                     .modelIndex=c(aaERactERERposERRARpos,aaERactERERposERRARneg),
                     .removeNonFunctionalRegs=F,
                     .nonFuncThreshold=0.3,
                     .exportplot=T,
                     .dendrogram = "none",
                     .Colv=F,
                     .Rowv = F,
                     setmaxto1RestZero = T,
                     .filename = "enhancerHeatmap.png",
                     .ColSideColors = aaacolside,
                     .RowSideCol = aarowsidecol
                     ,enhancerIndex = c(which(ER.associated.genes.MotifBased.reg.elements.Interest == 1), which(ER.associated.genes.MotifBased.reg.elements.Interest == 2))
                     )

ParameterHeatMap(parameterMat = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[c(aaERactERERposERRARpos,aaERactERERposERRARneg), ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row",
                 .RowSideColors = aaacolside)

PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output,
                    modelIndex = c(aaERactERERposERRARpos,aaERactERERposERRARneg),
                    .Colv = T,
                    .Rowv = T,
                    .dendrogram = "both",
                    .ColSideColors = aaacolside)

plot_ensemble_output(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output, ClusterNumber = 1, OneGeneManyModel = T,groupvsgroup = T, group1 = aaERactERERposERRARpos, group2 = aaERactERERposERRARneg, Rank = 2)

aa1 <- which.min(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output$averageSSE[aaERactERERposERRARpos])
aaERactERERposERRARpos[aa1]

aa2 <- which.min(Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output$averageSSE[aaERactERERposERRARneg])
aaERactERERposERRARneg[aa2]

CompareModelPlot(GEMSTAT_output_Reader_ensemble.Output = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output,
                 LogFileReader.Output = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param,
                 EnhWeightReader.Output = Exp18.GEMSTAT.All.Filtered.sigmoid_no_objss_5$enh,
                 modelIndex = c(157,159),
                 ..col_vec = col_vector,
                 .geneNames=rownames(ER.associated.genes.MotifBased.ExpMat),
                 ..exportplot = T,
                 plotDir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Ensemble/EXP18",
                 .setmaxto1RestZero = T)

########################################################################################################################
########################################################################################################################
######################   EXP 19_using -ln(sum(exp(- 5 * SSE(enhancer(j))))) as the obj functionr###############################################################
######################   using genes found in motif analysis
######################   Adding JUN motif (only JUN_2)[will add JUN_1 in the next experiment]
ER.associated.genes.MotifAnalysis.Timeseries.Exp.Mat
aa1 <- match(ER.associated.genes.ERERgt2.ERRAReq0.GEMSTAT, rownames(GSE78167.RNAseq.Avg.Norm))
aa2 <- match(ER.associated.genes.ERERgt3.ERRARgt3.GEMSTAT, rownames(GSE78167.RNAseq.Avg.Norm))
aa1 <- aa1[!is.na(aa1)]
aa2 <- aa2[!is.na(aa2)]

sum(ER.associated.genes.ERERgt2.ERRAReq0.GEMSTAT %in% rownames(GSE78167.RNAseq.Diff))
sum(ER.associated.genes.ERERgt3.ERRARgt3.GEMSTAT %in% rownames(GSE78167.RNAseq.Diff))

GSE78167.RNAseq.Diff

par(mfrow = c(5, 5), mar= c(0.1, 0.1, 0.1, 0.1))
for(i in 1:length(aa1)){
  plot(GSE78167.RNAseq.Avg.Norm[aa1[i], ], type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main="")
}
par(mfrow = c(5, 5), mar= c(0.1, 0.1, 0.1, 0.1))
for(i in 1:length(aa2)){
  plot(GSE78167.RNAseq.Avg.Norm[aa2[i], ], type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main="")
}
ER.associated.genes.MotifBased <- c(ER.associated.genes.ERERgt2.ERRAReq0.GEMSTAT, ER.associated.genes.ERERgt3.ERRARgt3.GEMSTAT)
#creating the expression matrix:
aa1 <- match(ER.associated.genes.MotifBased, rownames(GSE78167.RNAseq.Avg.Norm))
ER.associated.genes.MotifBased.ExpMat <- GSE78167.RNAseq.Avg.Norm[aa1, ]

#creating the motifs:
TF.motifs.TimeSeries.loose.plJUN2 <- lapply(X = c(TF.motifs.TimeSeries.loose.t, TF.motifs.TimeSeries.AP1.t[2]), t) 
TF.motifs.TimeSeries.loose.plJUN2.count <- lapply(TF.motifs.TimeSeries.loose.plJUN2, PWMtoCount)

#creating TF_expression
#JUN ENTREZ ID: "3725"
TF.Expression.TimeSeries

aa2 <- match("3725", rownames(GSE78167.RNAseq.Avg.Norm))
TF.Expression.TimeSeries.plJUN_2 <- rbind(TF.Expression.TimeSeries, GSE78167.RNAseq.Avg.Norm[aa2, ])
rownames(TF.Expression.TimeSeries.plJUN_2)[3] <- "JUN_2"

#adding basal expression to the parameter sampling space
#change the range of parameters in GEMSTAT:
# double ExprPar::min_weight = 0.01;
# double ExprPar::max_weight = 1000;
# double ExprPar::min_interaction = 0.01;	
# double ExprPar::max_interaction = 1000;
# double ExprPar::min_effect_Thermo = 0.01;	
# double ExprPar::max_effect_Thermo = 1000;
# double ExprPar::min_repression = 1.0E-5;
# double ExprPar::max_repression = 1; 
# double ExprPar::min_basal_Thermo = 1.0E-4;	
# double ExprPar::max_basal_Thermo = 1;
#creating range matrix
aa <- matrix(0L, nrow = 9, ncol = 3)
aa[1,] <- c(-2, 1, 1) #ESR1_Binding
aa[2,] <- c(0,  3, 1) #ESR1_txp
aa[3,] <- c(-2, 1, 1) #RARA_binding
aa[4,] <- c(-2, 3, 1) #RARA_txp
aa[5,] <- c(-2, 3, 1) #JUN_2_binding
aa[6,] <- c(-2, 3, 1) #JUN_2_txp
aa[7,] <- c(-4, 0, 1) #basal expression
aa[8,] <- c(1,  3, 1) #ESR1_ESR1
aa[9,] <- c(-2, 3, 1) #RARA_ESR1

#ex
# creating input for GEMSTAT
GEMSTAT_input_constructor_Ensemble(clusterIndex = list(c(1:nrow(ER.associated.genes.MotifBased.ExpMat))),
                                   gene_expression_Mat = ER.associated.genes.MotifBased.ExpMat, 
                                   .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                   .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                   .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                   .promoter.sequence=Promoter.Sequence.Char,
                                   .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                   .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_Max/seq2expr",
                                   .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                   experiment.nu = 19,
                                   #.par.exp.nu = 17,
                                   .motifs=TF.motifs.TimeSeries.loose.plJUN2.count,
                                   .TF_expression=TF.Expression.TimeSeries.plJUN_2,
                                   .max.enh.length = numeric(0),
                                   ..o="DIRECT", coopertingTFIndex=rbind(c(1,1),c(1,2))
                                   , ..oo="SSE", ..ct=11, ..sigma=character(0),
                                   ..na=10
                                   ,.par.RangeMat = aa,
                                   .nu.sample=20,
                                   writeParameters = T)
#reading output
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5 <- list()
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]] <- GEMSTAT_output_Reader_ensemble(.exp.nu = 19, .number.Of.Clusters = 1, copyOutput = F)
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]] <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_19/Cluster_1/Outputs",
                                                                    removeDups =F, nu.TFs = 3, nu.coop.pairs = 2,
                                                                    parameterNames = c("ER_bind", "ER_txp", "RAR_bind", "RAR_txp", "JUN2_bind", "JUN2_txp", "basal","ER_ER","ER_RAR"))
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[3]] <- EnhWeightReader("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_19/Cluster_1/Outputs")
names(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5) <- c("output", "param", "enh")

aad <- sort(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE, decreasing = F, index.return = T)$ix 

plot_ensemble_output(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                     ClusterNumber = 1,
                     paramnumber = 326)
plot_ensemble_output(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                     ClusterNumber = 1,
                     paramnumber = 174)
plot_ensemble_output(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                     ClusterNumber = 1,
                     OneGeneManyModel = T, Rank = 1)
plot_ensemble_output(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                     ClusterNumber = 1,
                     OneGeneManyModel = T, GeneNumber = 61)
plot_ensemble_output(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                     ClusterNumber = 1,
                     OneGeneManyModel = T,
                     plotAmong = c(150, 38, 45, 93, 39), varAmong = c(150, 38, 45, 93, 39), Rank = 3)

######################Filter the similar models
aaa <- logSimilarityCheck(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]], logscale = T, thresh = 0.1)
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_filterpar <- aaa
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5 <- FilterForModels(ModelIndex = setdiff(c(1:length(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE)), Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_filterpar$removed),
                                                                 GEMSTAT_output_Reader_ensemble.output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                                                                 EnhWeightReader.output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[3]],
                                                                 basedOnSimilarityCheck = F,
                                                                 .LogFileReader.Output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]])
remove(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_filterpar)
aad <- sort(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE, decreasing = F, index.return = T)$ix 

par(mfrow = c(1,1), mar = c(4,4,4,4))
plot(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE[aad], main = "obj_ss_5_performance_Filtered", ylab = "SSE", xlab="models")
abline(h = seq(0.160,0.230, 0.001), col = "grey", lty = 3, lwd = 0.7)
#########
aaarowsidecol <- c(rep("blue", 22), rep("pink", 22))
PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                    modelIndex = aad,
                    .Colv = F,
                    .Rowv = F,
                    .dendrogram = "none",
                    .RowSideColors = aaarowsidecol)
aadd <- which(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE < 0.18)
ParameterHeatMap(parameterMat = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[aadd, ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row")
###################################################
#create the character vector of the sequence of all reg elements of the modeled genes:
aaa <- match(rownames(ER.associated.genes.MotifBased.ExpMat), names(Vicinity100kb.Enhancer.plusPromoter.By.gene))
ER.associated.genes.MotifBased.reg.elements.seqANDlength <- WriteFastaOfBag(enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene[aaa],
                                                                            Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                                                            promoter.promoter.int=Promoter.gene.PromoterIntList[aaa],
                                                                            promoter.sequence=Promoter.Sequence.Char,
                                                                            return.Length = T,
                                                                            returnListNotwrite = T)
ER.associated.genes.MotifBased.reg.elements.seq <- do.call(c, ER.associated.genes.MotifBased.reg.elements.seq$sequence)

##sequence of regulatory elements of interest:
ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence
ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence
# ER.associated.genes.MotifBased.reg.elements.Interest is an integer vector with one entry associated with each reg element. entry is 1: ER no RAR, entry is 2: ER and RAR, otherwise zero
ER.associated.genes.MotifBased.reg.elements.Interest <- integer(length = length(ER.associated.genes.MotifBased.reg.elements.seq))
#Find which of the reg elements are the ones of interest, assign color blue to ones without RAR and pink to ones with RAR, and white to others
aarowsidecol <- rep("white", length(ER.associated.genes.MotifBased.reg.elements.seq))
aa <- which(toupper(ER.associated.genes.MotifBased.reg.elements.seq) %in% toupper(ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence))
aarowsidecol[aa] <- "pink"
ER.associated.genes.MotifBased.reg.elements.Interest[aa] <- 2
aa <- which(toupper(ER.associated.genes.MotifBased.reg.elements.seq) %in% toupper(ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence))
aarowsidecol[aa] <- "blue"
ER.associated.genes.MotifBased.reg.elements.Interest[aa] <- 1
###################################################################
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[3]],
                     .removeNonFunctionalRegs = T, .nonFuncThreshold = 0.1,
                     .exportplot = T,
                     .modelIndex = aad,
                     .dendrogram = "column",
                     .Colv = T, .Rowv = F,
                     setmaxto1RestZero = T
                     , .RowSideCol = aarowsidecol
)
GSE78167.RNAseq.Avg.ChosenGenes.Norm.Reg.MotifBased.Coor <- RegulatoryCoordinateExtractor(gene.Names = rownames(ER.associated.genes.MotifBased.ExpMat))
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_RegElemWeightCoordinate <- RegElemWeightCoordinate(RegulatoryCoordinateExtractorOutput = GSE78167.RNAseq.Avg.ChosenGenes.Norm.Reg.MotifBased.Coor,
                                                                                                 EnhWeightReaderOutput = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[3]],
                                                                                                 plotheatmap = F,
                                                                                                 setmaxto1RestZero = T)
####Look which models use the enhancers of interest the most: THink how to plot this??
# for each gene compare the average performance of models which use enhancer of interest vs other enhancers
aa <- integer(length(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[3]]))
for( i in 1: length(aa)){
  aa[i] <- sum(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_RegElemWeightCoordinate[[5 + i]][which(ER.associated.genes.MotifBased.reg.elements.Interest > 0)])
}
hist(aa ,main = "no. of enhancers of interest used/out of 61")

###Percent gene explained
aaPR <- PredictedvsRealCorRMSD(GEMSTAT_output_Reader_ensemble.output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                               corrMethod = "pearson" )

aaa <-      PlotPercentGenes(..GEMSTAT_output_Reader_ensemble.output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]],
                             .LogFileReader.output=Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]],
                             .EnhWeightReader.output=Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[3]],
                             RMSDrange = c(0,0.15,0.2,1),
                             CorRange = c(1, 0.8, 0.6, -1),
                             .corrMethod = "pearson",
                             Num.Cluster = 4,
                             sortBasedOn = 1,
                             .CorRMSDList = aaPR,
                             exportplot = T)

#find models with ER as activator, ER-ER interaction positive
aaERact <- which(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[, 2] > 1)
aaERERpos <- which(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[, 8] > 1)
aaERRARpos <- which(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[, 9] > 1)
aaERRARneg <- which(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[, 9] < 1)
aaJunInh <- which(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[, 6] < 1)
aaERactERERpos <- intersect(aaERact, aaERERpos)
aaERactERERposJunInh <- intersect(aaERactERERpos, aaJunInh)

aaERactERERposERRARpos <- intersect(aaERactERERpos, aaERRARpos)
aaERactERERposERRARneg <- intersect(aaERactERERpos, aaERRARneg)
hist(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE[aaERactERERpos])
aadd <- which(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[1]]$averageSSE < 0.18)

ParameterHeatMap(parameterMat = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[aaERactERERpos, ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row")

ParameterHeatMap(parameterMat = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5[[2]]$parameters[aaERactERERposJunInh, ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row")

aaacolside <- c(rep("pink", length(aaERactERERposERRARpos)), rep("blue", length(aaERactERERposERRARneg)))
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5$enh,
                     .col_vec=col_vector,
                     .modelIndex=c(aaERactERERposERRARpos,aaERactERERposERRARneg),
                     .removeNonFunctionalRegs=F,
                     .nonFuncThreshold=0.3,
                     .exportplot=T,
                     .dendrogram = "none",
                     .Colv=F,
                     .Rowv = F,
                     setmaxto1RestZero = T,
                     .filename = "enhancerHeatmap.png",
                     .ColSideColors = aaacolside,
                     .RowSideCol = aarowsidecol
                     ,enhancerIndex = c(which(ER.associated.genes.MotifBased.reg.elements.Interest == 1), which(ER.associated.genes.MotifBased.reg.elements.Interest == 2))
)

ParameterHeatMap(parameterMat = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param$parameters[c(aaERactERERposERRARpos,aaERactERERposERRARneg), ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row",
                 .RowSideColors = aaacolside)

PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output,
                    modelIndex = c(aaERactERERposERRARpos,aaERactERERposERRARneg),
                    .Colv = T,
                    .Rowv = T,
                    .dendrogram = "both",
                    .ColSideColors = aaacolside)

plot_ensemble_output(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output, ClusterNumber = 1, OneGeneManyModel = T,groupvsgroup = T, group1 = aaERactERERposERRARpos, group2 = aaERactERERposERRARneg, Rank = 2)

aa1 <- which.min(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output$averageSSE[aaERactERERposERRARpos])
aaERactERERposERRARpos[aa1]

aa2 <- which.min(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output$averageSSE[aaERactERERposERRARneg])
aaERactERERposERRARneg[aa2]

CompareModelPlot(GEMSTAT_output_Reader_ensemble.Output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5$output,
                 LogFileReader.Output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5$param,
                 EnhWeightReader.Output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5$enh,
                 modelIndex = c(157,159),
                 ..col_vec = col_vector,
                 .geneNames=rownames(ER.associated.genes.MotifBased.ExpMat),
                 ..exportplot = T,
                 plotDir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Ensemble/EXP18",
                 .setmaxto1RestZero = T)
########################################################################################################################
########################################################################################################################
######################   EXP 19_using -ln(sum(exp(- 5 * SSE(enhancer(j))))) as the obj functionr###############################################################
######################   using genes found in motif analysis
######################   Adding JUN motif (only JUN_2)[will add JUN_1 in the next experiment]
######################   Using logSumExp trick

#reading output
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick <- list()
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]] <- GEMSTAT_output_Reader_ensemble(.exp.nu = 19, .number.Of.Clusters = 1, copyOutput = F)
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[2]] <- LogFileReader(logdir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_19/Cluster_1/Outputs",
                                                                    removeDups =F, nu.TFs = 3, nu.coop.pairs = 2,
                                                                    parameterNames = c("ER_bind", "ER_txp", "RAR_bind", "RAR_txp", "JUN2_bind", "JUN2_txp", "basal","ER_ER","ER_RAR"))
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[3]] <- EnhWeightReader("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble/Experiment_19/Cluster_1/Outputs")
names(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick) <- c("output", "param", "enh")

aad <- sort(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]]$averageSSE, decreasing = F, index.return = T)$ix 

plot_ensemble_output(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]],
                     ClusterNumber = 1,
                     paramnumber = aad[1])
plot_ensemble_output(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]],
                     ClusterNumber = 1,
                     paramnumber = 174)
plot_ensemble_output(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]],
                     ClusterNumber = 1,
                     OneGeneManyModel = T, Rank = 1)
plot_ensemble_output(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]],
                     ClusterNumber = 1,
                     OneGeneManyModel = T, GeneNumber = 61)
plot_ensemble_output(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]],
                     ClusterNumber = 1,
                     OneGeneManyModel = T,
                     plotAmong = c(150, 38, 45, 93, 39), varAmong = c(150, 38, 45, 93, 39), Rank = 3)

######################Filter the similar models
aaa <- logSimilarityCheck(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[2]], logscale = T, thresh = 0.1)
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_filterpar_trick <- aaa
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick <- FilterForModels(ModelIndex = setdiff(c(1:length(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]]$averageSSE)), Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_filterpar_trick$removed),
                                                                 GEMSTAT_output_Reader_ensemble.output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]],
                                                                 EnhWeightReader.output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[3]],
                                                                 basedOnSimilarityCheck = F,
                                                                 .LogFileReader.Output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[2]])
remove(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_filterpar_trick)
aad <- sort(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]]$averageSSE, decreasing = F, index.return = T)$ix 

par(mfrow = c(1,1), mar = c(4,4,4,4))
plot(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]]$averageSSE[aad], main = "obj_ss_5_performance_Filtered", ylab = "SSE", xlab="models")
abline(h = seq(0.160,0.230, 0.001), col = "grey", lty = 3, lwd = 0.7)
###############
aaarowsidecol <- c(rep("blue", 22), rep("pink", 22))
PerformanceHeattMap(GEMSTAT_output_Reader_ensemble.Output =Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]],
                    modelIndex = aad,
                    .Colv = F,
                    .Rowv = F,
                    .dendrogram = "none",
                    .RowSideColors = aaarowsidecol)
aadd <- which(Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]]$averageSSE < 0.18)
ParameterHeatMap(parameterMat = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[2]]$parameters[aadd, ],
                 exportplot = T,
                 .Rowv = T,
                 .dendrogram = "row")
aarowsidecol <- rep("white", length(ER.associated.genes.MotifBased.reg.elements.seq))
aa <- which(toupper(ER.associated.genes.MotifBased.reg.elements.seq) %in% toupper(ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence))
aarowsidecol[aa] <- "pink"
aa <- which(toupper(ER.associated.genes.MotifBased.reg.elements.seq) %in% toupper(ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence))
aarowsidecol[aa] <- "blue"
###################################################################
RegElemWeightHeatMap(.EnhWeightReaderOutput = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[3]],
                     .removeNonFunctionalRegs = T, .nonFuncThreshold = 0.1,
                     .exportplot = T,
                     .modelIndex = aad,
                     .dendrogram = "column",
                     .Colv = T, .Rowv = F,
                     setmaxto1RestZero = T
                     , .RowSideCol = aarowsidecol
)
GSE78167.RNAseq.Avg.ChosenGenes.Norm.Reg.MotifBased.Coor <- RegulatoryCoordinateExtractor(gene.Names = rownames(ER.associated.genes.MotifBased.ExpMat))
Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_RegElemWeightCoordinate_trick <- RegElemWeightCoordinate(RegulatoryCoordinateExtractorOutput = GSE78167.RNAseq.Avg.ChosenGenes.Norm.Reg.MotifBased.Coor,
                                                                                                 EnhWeightReaderOutput = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[3]],
                                                                                                 plotheatmap = F,
                                                                                                 setmaxto1RestZero = T)
###Percent gene explained
aaPR <- PredictedvsRealCorRMSD(GEMSTAT_output_Reader_ensemble.output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]],
                               corrMethod = "pearson" )

aaa <-      PlotPercentGenes(..GEMSTAT_output_Reader_ensemble.output = Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[1]],
                             .LogFileReader.output=Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[2]],
                             .EnhWeightReader.output=Exp19.GEMSTAT.All.Filtered.sigmoid_no_objss_5_trick[[3]],
                             RMSDrange = c(0,0.15,0.2,1),
                             CorRange = c(1, 0.8, 0.6, -1),
                             .corrMethod = "pearson",
                             Num.Cluster = 4,
                             sortBasedOn = 1,
                             .CorRMSDList = aaPR,
                             exportplot = T)

########################################################################################################################
########################################################################################################################
######################   EXP 20_using rmse/cor as the obj function#######################
######################   using genes found in motif analysis
######################   using only one enhancer per gene --> not multi enhancer
######################   Adding JUN motif (only JUN_2)

#given the enhacners associated with the two categories: ERonly and ER-RAR, for each enhacner find one gene associated with the enhancer
ER.associated.reg.elements.ERERRAR.Filtered.sequence.LT10KB$ERgt1RARgt1
ER.associated.reg.elements.ERERRAR.Filtered.sequence.LT10KB$ERgt2RAReq0
# full sequences and length of regulatory elements for the 44 genes of interest: first 22 only ER, second 22 both
ER.associated.genes.MotifBased.reg.elements.seqANDlength

#create an charachter vector with one reg element per gene: that reg element should be in ER.associated.reg.elements.ERERRAR.Filtered.sequence.LT10KB
ER.associated.genes.MotifBased.reg.elements.SingleEntry <- character(length(ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence))
for(i in 1:22){
  aadone=0
  for(j in 1:length(ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence[[i]])){
    if(ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence[[i]][j] %in% ER.associated.reg.elements.ERERRAR.Filtered.sequence.LT10KB$ERgt2RAReq0){
      ER.associated.genes.MotifBased.reg.elements.SingleEntry[i] <- ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence[[i]][j]
      aadone=1
      break
    }
  }
  if (aadone==0){
    for(j in 1:length(ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence[[i]])){
      if(ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence[[i]][j] %in% ER.associated.reg.elements.ERERRAR.Filtered.sequence$ERgt2RAReq0){
        ER.associated.genes.MotifBased.reg.elements.SingleEntry[i] <- ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence[[i]][j]
        aadone=1
        break
      }
    }
  }
}

for(i in 23:44){
  aadone=0
  for(j in 1:length(ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence[[i]])){
    if(ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence[[i]][j] %in% ER.associated.reg.elements.ERERRAR.Filtered.sequence.LT10KB$ERgt1RARgt1){
      ER.associated.genes.MotifBased.reg.elements.SingleEntry[i] <- ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence[[i]][j]
      aadone=1
      break
    }
  }
  if (aadone==0){
    for(j in 1:length(ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence[[i]])){
      if(ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence[[i]][j] %in% ER.associated.reg.elements.ERERRAR.Filtered.sequence$ERgt1RARgt1){
        ER.associated.genes.MotifBased.reg.elements.SingleEntry[i] <- ER.associated.genes.MotifBased.reg.elements.seqANDlength$sequence[[i]][j]
        aadone=1
        break
      }
    }
  }
}
########################################################################################################################
ER.associated.genes.MotifBased.reg.elements.SingleEntry # is the single enhacner per gene character vector
########################################################################################################################
GEMSTAT_input_constructor_Ensemble(clusterIndex = list(c(1:nrow(ER.associated.genes.MotifBased.ExpMat))),
                                   gene_expression_Mat = ER.associated.genes.MotifBased.ExpMat, 
                                   .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                   .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                   .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                   .promoter.sequence=Promoter.Sequence.Char,
                                   singleEnhancer = T,
                                   Manual_enhancer_Seq = ER.associated.genes.MotifBased.reg.elements.SingleEntry,
                                   .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                   .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_Max/seq2expr",
                                   .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                   experiment.nu = 20,
                                   .par.exp.nu = 19,
                                   .motifs=TF.motifs.TimeSeries.loose.plJUN2.count,
                                   .TF_expression=TF.Expression.TimeSeries.plJUN_2,
                                   .max.enh.length = numeric(0),
                                   ..o="DIRECT", coopertingTFIndex=rbind(c(1,1),c(1,2))
                                   , ..oo="SSE", ..ct=11, ..sigma=character(0),
                                   ..na=10
                                   ,.par.RangeMat = aa,
                                   .nu.sample=20,
                                   writeParameters = F)

########################################################################################################################
########################################################################################################################
######################   EXP 21_using -ln(sum(exp(- 5 * SSE(enhancer(j))))) as the obj functionr#######################
######################   Using logSumExp trick for numerical stability
######################   using genes found in motif analysis
######################   Adding JUN motif (both JUN_1 and JUN_2)
#creating the motifs:
TF.motifs.TimeSeries.loose.plJUN.all <- lapply(X = c(TF.motifs.TimeSeries.loose.t, TF.motifs.TimeSeries.AP1.t), t) 
TF.motifs.TimeSeries.loose.plJUN.all.count <- lapply(TF.motifs.TimeSeries.loose.plJUN.all, PWMtoCount)

#creating TF_expression
#JUN ENTREZ ID: "3725"
TF.Expression.TimeSeries

aa2 <- match("3725", rownames(GSE78167.RNAseq.Avg.Norm))
TF.Expression.TimeSeries.plJUN_all <- rbind(TF.Expression.TimeSeries, GSE78167.RNAseq.Avg.Norm[aa2, ], GSE78167.RNAseq.Avg.Norm[aa2, ])

rownames(TF.Expression.TimeSeries.plJUN_all)[3:4] <- c("JUN_1", "JUN_2")

#adding basal expression to the parameter sampling space
#change the range of parameters in GEMSTAT:
# double ExprPar::min_weight = 0.01;
# double ExprPar::max_weight = 1000;
# double ExprPar::min_interaction = 0.01;	
# double ExprPar::max_interaction = 1000;
# double ExprPar::min_effect_Thermo = 0.01;	
# double ExprPar::max_effect_Thermo = 1000;
# double ExprPar::min_repression = 1.0E-5;
# double ExprPar::max_repression = 1; 
# double ExprPar::min_basal_Thermo = 1.0E-4;	
# double ExprPar::max_basal_Thermo = 1;
#creating range matrix

aa <- matrix(0L, nrow = 11, ncol = 3)
aa[1,] <- c(-2, 1, 1) #ESR1_Binding
aa[2,] <- c(0,  3, 1) #ESR1_txp
aa[3,] <- c(-2, 1, 1) #RARA_binding
aa[4,] <- c(-2, 3, 1) #RARA_txp
aa[5,] <- c(-2, 3, 1) #JUN_1_binding
aa[6,] <- c(-2, 3, 1) #JUN_1_txp
aa[7,] <- c(-2, 3, 1) #JUN_2_binding
aa[8,] <- c(-2, 3, 1) #JUN_2_txp
aa[9,] <- c(-4, 0, 1) #basal expression
aa[10,] <- c(1,  3, 1) #ESR1_ESR1
aa[11,] <- c(-2, 3, 1) #RARA_ESR1

GEMSTAT_input_constructor_Ensemble(clusterIndex = list(c(1:nrow(ER.associated.genes.MotifBased.ExpMat))),
                                   gene_expression_Mat = ER.associated.genes.MotifBased.ExpMat, 
                                   .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                   .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                   .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                   .promoter.sequence=Promoter.Sequence.Char,
                                   .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Ensemble",
                                   .GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_Max/seq2expr",
                                   .sharedmountsdir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/Ensemble/",
                                   experiment.nu = 20,
                                   #.par.exp.nu = 17,
                                   .motifs=TF.motifs.TimeSeries.loose.plJUN.all.count,
                                   .TF_expression=TF.Expression.TimeSeries.plJUN_all,
                                   .max.enh.length = numeric(0),
                                   ..o="DIRECT", coopertingTFIndex=rbind(c(1,1),c(1,2))
                                   , ..oo="SSE", ..ct=11, ..sigma=character(0),
                                   ..na=10
                                   ,.par.RangeMat = aa,
                                   .nu.sample=20,
                                   writeParameters = T)


