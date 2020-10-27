#experiments:
####################
#I want to organize the kd and treatment datasets into several experiments each with one control experiment
colnames(GSE73663.Inh.RNAseq) #1 Zeynep     ctrl: 11,12 
#1"AZD334h1"   #treatment time: 4h inhibitor: 1 uM AZD6422 treatment: 1 uM PaPE
#2"AZD334h2"   #treatment time: 4h inhibitor: 1 uM AZD6422 treatment: 1 uM PaPE
#3"AZDE24h1"   #treatment time: 4h inhibitor: 1 uM AZD6422 treatment: 10 nM E2
#4"AZDE24h2"   #treatment time: 4h inhibitor: 1 uM AZD6422 treatment: 10 nM E2
#5"AZDVeh4h1"  #treatment time: 4h inhibitor: 1 uM AZD6422 treatment: Veh
#6"AZDVeh4h2"  #treatment time: 4h inhibitor: 1 uM AZD6422 treatment: Veh
#7"Ctrl334h1"  #treatment time: 4h inhibitor: veh          treatment: 1 uM PaPE
#8"Ctrl334h2"  #treatment time: 4h inhibitor: veh          treatment: 1 uM PaPE
#9"CtrlE24h1"  #treatment time: 4h inhibitor: veh          treatment: 10 nM E2
#10"CtrlE24h2"  #treatment time: 4h inhibitor: veh          treatment: 10 nM E2
#11"CtrlVeh4h1" #treatment time: 4h inhibitor: veh          treatment: Veh
#12"CtrlVeh4h2" #treatment time: 4h inhibitor: veh          treatment: Veh
#13"PP242334h1" #treatment time: 4h inhibitor: PP242        treatment: 1 uM PaPE
#14"PP242334h2" #treatment time: 4h inhibitor: PP242        treatment: 1 uM PaPE
#15"PP242E24h1" #treatment time: 4h inhibitor: PP242        treatment: 10 nM E2
#16"PP242E24h2" #treatment time: 4h inhibitor: PP242        treatment: 10 nM E2
#17"PP242Veh4h1"#treatment time: 4h inhibitor: PP242        treatment: Veh
#18"PP242Veh4h2"#treatment time: 4h inhibitor: PP242        treatment: Veh
#I will organize them such that the control experiment for the inhibitors is that inhibitor plus veh
#exp1 : control: 5"AZDVeh4h1" other : 1 "AZD334h1", 3 AZDE24h1
#exp2 : control: 6"AZDVeh4h2" other : 2"AZD334h2",4 AZDE24h2
#exp3 : control: 17 PP242Veh4h1 other: 13 PP242334h1 , 15 PP242E24h1
#exp4 : control: 18 PP242Veh4h2 other: 14 PP242334h2 ,16 PP242E24h2
#exp5 : control: 11 CtrlVeh4h1 other:  7 Ctrl334h1 , 9 CtrlE24h1
#exp6 : control: 12 CtrlVeh4h2 other:  8 Ctrl334h2 , 10 CtrlE24h2
aa <- c("MAPK_Veh_1", "MAPK_PaPE_1", "MAPK_E2_1",
        "MAPK_Veh_2", "MAPK_PaPE_2", "MAPK_E2_2",
        "mTOR_Veh_1", "mTOR_PaPE_1", "mTOR_E2_1",
        "mTOR_Veh_2", "mTOR_PaPE_2", "mTOR_E2_2",
        "Veh_Veh_1", "Veh_PaPE_1", "Veh_E2_1",
        "Veh_Veh_2", "Veh_PaPE_2", "Veh_E2_2")
####################################################################################################
####################################################################################################
####################################################################################################
colnames(GSE86316.RNAseq) #2 Men1 kd   ctrl: 1,2
"kd_smpl_1_s1_treat"               #Veh_ShCtrl_1
"kd_smpl_2_s1_treat"               #Veh_ShCtrl_2
"MCF7_E2_rep1_s1_treat"            #E2_ShCtrl_1 10nM estradiol 45 mins
"MCF7_E2_rep2_s1_treat"            #E2_ShCtrl_2 10nM estradiol 45 mins
"MCF7_Veh_KD_MEN1_rep1_s1_treat"   #Veh_ShMEN1_1
"MCF7_Veh_KD_MEN1_rep2_s1_treat"   #Veh_ShMEN1_2
"MCF7_E2_KD_MEN1_rep1_s1_treat"    #E2_ShMEN1_1 10nM estradiol 45 mins
"kd_smpl_8_s1_treat"               #E2_ShMEN1_2 10nM estradiol 45 mins
####################
#exp1 : control: Veh_ShCtrl_1        other: E2_ShCtrl_1 
#exp2 : control: Veh_ShCtrl_2        other: E2_ShCtrl_2
#exp3 : control: Veh_ShMEN1_1        other: E2_ShMEN1_1
#exp4 : control: Veh_ShMEN1_2        other: E2_ShMEN1_2
####################################################################################################
####################################################################################################
####################################################################################################
#Drop CTCF for now. Since the changes causes by it knockdown are hard to model
colnames(mutualGSE53532RNA) ### CTCF kd     ctrl: 1 #the name of dataset is wrong, the correct name is GSM2257523
#"siNT_E0h"   
#"siNT_E3h"   
#"siCTCF_E0h" 
#"siCTCF_E3h"
#############
#exp1: control: siNT_E0h      other: siNT_E3h
#exp2: control: siCTCF_E0h    other: siCTCF_E3h
####################################################################################################
####################################################################################################
####################################################################################################
colnames(GSE55922.RNAseq) #3 BRD4 kd     cntrl: 5,6
#"brd4.e2a"     #siBRD4 -- Vehicle
#"brd4.e2b"     #siBRD4 -- Vehicle
#"brd4pluse2a"  #siBRD4 -- 10 nmol/L 17-β-estradiol (Sigma-Aldrich) for 2 hours.
#"brd4pluse2b"  #siBRD4 -- 10 nmol/L 17-β-estradiol (Sigma-Aldrich) for 2 hours.
#"cont.e2a"     #sictrl -- Vehicle
#"cont.e2b"     #sictrl -- Vehicle
#"contpluse2a"  #sictrl -- 10 nmol/L 17-β-estradiol (Sigma-Aldrich) for 2 hours.
#"contpluse2b"  #sictrl -- 10 nmol/L 17-β-estradiol (Sigma-Aldrich) for 2 hours.
#"jq1.e2a"      #treated with JQ1 which is a BRD4 inhibitor -- Vehicle
#"jq1.e2b"      #treated with JQ1 which is a BRD4 inhibitor -- Vehicle
#"jq1pluse2a"   #treated with JQ1 which is a BRD4 inhibitor -- 10 nmol/L 17-β-estradiol (Sigma-Aldrich) for 2 hours.
#"jq1pluse2b"  #treated with JQ1 which is a BRD4 inhibitor -- 10 nmol/L 17-β-estradiol (Sigma-Aldrich) for 2 hours.
##############
#exp1: control: brd4.e2a other: brd4pluse2a
#exp2: control: brd4.e2b other: brd4pluse2b
#exp3: control: cont.e2a other: contpluse2a
#exp4: control: cont.e2b other: contpluse2b
#exp5: control: jq1.e2a other: jq1pluse2a
#exp6: control: jq1.e2b other: jq1pluse2b
aa <- c("siBRD4_Veh_1", "siBRD4_E2_1",
        "siBRD4_Veh_2", "siBRD4_E2_2",
        "sictrl_Veh_1", "sictrl_E2_1",
        "sictrl_Veh_2", "sictrl_E2_2",
        "JQ1_Veh_1", "JQ1_E2_1",
        "JQ1_Veh_2", "JQ1_E2_2")

####################################################################################################
####################################################################################################
####################################################################################################
colnames(GSE67295.RNAseq) #4 #E2 and other treatments ctrl : 11,12
# [1] "Sample19.RNASeq.Veh.Expt1.Rep1.fastq.gz"  1      treatment: Veh
# [2] "Sample20.RNASeq.Veh.Expt1.Rep2.fastq.gz"  1      treatment: Veh
# [3] "Sample21.RNASeq.E2.Rep1.fastq.gz"         1      treatment: E2 3hours
# [4] "Sample22.RNASeq.E2.Rep2.fastq.gz"         1      treatment: E2 3hours
# [5] "Sample23.RNASeq.E2.ICI.Rep1.fastq.gz"     1      treatment: Fulvestrant + E2 3hours
# [6] "Sample24.RNASeq.E2.ICI.Rep2.fastq.gz"     1      treatment: Fulvestrant + E2 3hours
# [7] "Sample25.RNASeq.IL1b.Rep1.fastq.gz"       1      treatment: Il1b 3hours
# [8] "Sample26.RNASeq.IL1b.Rep2.fastq.gz"       1      treatment: Il1b 3hours
# [9] "Sample27.RNASeq.IL1b.ICI.Rep1.fastq.gz"   1      treatment: Il1b + Fulvestrant 3hrs
# [10] "Sample28.RNASeq.IL1b.ICI.Rep2.fastq.gz"  1      treatment: Il1b + Fulvestrant 3hrs
# [11] "Sample29.RNASeq.Veh.Rep1.fastq.gz"       1      treatment: Veh
# [12] "Sample30.RNASeq.Veh.Rep2.fastq.gz"       1      treatment: Veh
# [13] "Sample31.RNASeq.TNF.Rep1.fastq.gz"       1      treatment: TNFa for 3h
# [14] "Sample32.RNASeq.TNF.Rep2.fastq.gz"       1      treatment: TNFa for 3h
# [15] "Sample33.RNASeq.TNF.ICI.Rep1.fastq.gz"         treatment: TNFa + Fulvestrant for 3h
# [16] "Sample34.RNASeq.TNF.ICI.Rep2.fastq.gz"         treatment: TNFa + Fulvestrant for 3h
# [17] "Sample35.RNASeq.Veh.Expt3.Rep1.fastq.gz"       treatment: Veh
# [18] "Sample36.RNASeq.Veh.Expt3.Rep2.fastq.gz"       treatment: Veh
# [19] "Sample37.RNASeq.E2.Expt3.Rep1.fastq.gz"        treatment: E2 3hours
# [20] "Sample38.RNASeq.E2.Expt3.Rep2.fastq.gz"        treatment: E2 3hours
# [21] "Sample39.RNASeqE2.TOT.Expt3.Rep1.fastq.gz"     treatment: 4OH-tamoxifen +E2 3hours 
# [22] "Sample40.RNASeq.E2.TOT.Expt3.Rep2.fastq.gz"    treatment: 4OH-tamoxifen +E2 3hours 
# [23] "Sample41.RNASeq.E2.TOT.IL1b.Expt3.Rep1.fastq.gz"    treatment: 4OH-tamoxifen +E2 + IL1b 3hours 
# [24] "Sample42.RNASeq.E2.TOT.IL1b.Expt3.Rep2.fastq.gz"    treatment: 4OH-tamoxifen +E2 + IL1b 3hours 
# [25] "Sample43.RNASeq.E2.TOT.TNF.Expt3.Rep1.fastq.gz"     treatment: 4OH-tamoxifen +E2 + TNFa 3hours 
# [26] "Sample44.RNASeq.E2.TOT.TNF.Expt3.Rep2.fastq.gz"     treatment: 4OH-tamoxifen +E2 + TNFa 3hours 
############################
#exp1: control: "Sample19.RNASeq.Veh.Expt1.Rep1.fastq.gz" other Sample21.RNASeq.E2.Rep1.fastq.gz, Sample23.RNASeq.E2.ICI.Rep1.fastq.gz, Sample25.RNASeq.IL1b.Rep1.fastq.gz, Sample27.RNASeq.IL1b.ICI.Rep1.fastq.gz
#exp2: control: Sample20.RNASeq.Veh.Expt1.Rep2.fastq.gz   other Sample22.RNASeq.E2.Rep2.fastq.gz, Sample24.RNASeq.E2.ICI.Rep2.fastq.gz, Sample26.RNASeq.IL1b.Rep2.fastq.gz, Sample28.RNASeq.IL1b.ICI.Rep2.fastq.gz
#exp3: control: Sample29.RNASeq.Veh.Rep1.fastq.gz         other Sample31.RNASeq.TNF.Rep1.fastq.gz, Sample33.RNASeq.TNF.ICI.Rep1.fastq.gz
#exp4: control: Sample30.RNASeq.Veh.Rep2.fastq.gz         other Sample32.RNASeq.TNF.Rep2.fastq.gz, Sample34.RNASeq.TNF.ICI.Rep2.fastq.gz
#exp5: control: Sample35.RNASeq.Veh.Expt3.Rep1.fastq.gz   other Sample39.RNASeqE2.TOT.Expt3.Rep1.fastq.gz, Sample41.RNASeq.E2.TOT.IL1b.Expt3.Rep1.fastq.gz, Sample43.RNASeq.E2.TOT.TNF.Expt3.Rep1.fastq.gz
#exp6: control: Sample36.RNASeq.Veh.Expt3.Rep2.fastq.gz   other Sample40.RNASeq.E2.TOT.Expt3.Rep2.fastq.gz,Sample42.RNASeq.E2.TOT.IL1b.Expt3.Rep2.fastq.gz, Sample44.RNASeq.E2.TOT.TNF.Expt3.Rep2.fastq.gz
c("Veh_1_1", "E2_1_1", "Fulv+E2_1_1",  "Il1b_1_1", "Fulv+Il1b_1_1",
  "Veh_1_2", "E2_1_2", "Fulv+E2_1_2",  "Il1b_1_2", "Fulv+Il1b_1_2",
  "Veh_2_1", "TNFa_2_1", "Fulv+TNFa_2_1",
  "Veh_2_2", "TNFa_2_2", "Fulv+TNFa_2_2",
  "Veh_3_1",  "E2_3_1", "Tamxfn+E2_3_1", "Tamxfn+E2+IL1b_3_1", "Tamxfn+E2+TNFa_3_1",
  "Veh_3_2",  "E2_3_2", "Tamxfn+E2_3_2", "Tamxfn+E2+IL1b_3_2", "Tamxfn+E2+TNFa_3_2")
####################################################################################################
####################################################################################################
####################################################################################################

colnames(GSE68358.RNAseq)#5 #E2 or E2 plus something ctrl : 1,5,9,10,11,14,18,24
# [1] "GSM1669160_MCF7_Control_3.txt"   1 MCF7_E2_Rep3
# [2] "GSM1669161_MCF7_R5020_5.txt"      MCF7_E2+R5020_3hr_Rep5
# [3] "GSM1669162_MCF7_Prog_3.txt"      1 MCF7_E2+Progesterone_3hr_Rep3
# [4] "GSM1669163_MCF7_Prog_1.txt"      1 MCF7_E2+Progesterone_3hr_Rep1
# [5] "GSM1669164_MCF7_Control_5.txt"   1 MCF7_E2_Rep5
# [6] "GSM1669165_MCF7_Prog_4.txt"      1 MCF7_E2+Progesterone_3hr_Rep4
# [7] "GSM1669166_MCF7_R5020_4.txt"      MCF7_E2+R5020_3hr_Rep4
# [8] "GSM1669167_MCF7_R5020_1.txt"      MCF7_E2+R5020_3hr_Rep1
# [9] "GSM1669168_MCF7_Control_6.txt"   1 MCF7_E2_Rep6
# [10] "GSM1669169_MCF7_Control_8.txt"  1 MCF7_E2_Rep8
# [11] "GSM1669170_MCF7_Control_1.txt"  1 MCF7_E2_Rep1
# [12] "GSM1669171_MCF7_R5020_3.txt"     MCF7_E2+R5020_3hr_Rep3
# [13] "GSM1669172_MCF7_R5020_7.txt"     MCF7_E2+R5020_3hr_Rep7
# [14] "GSM1669173_MCF7_Control_7.txt"  1 MCF7_E2_Rep7
# [15] "GSM1669174_MCF7_Prog_2.txt"     1 MCF7_E2+Progesterone_3hr_Rep2
# [16] "GSM1669175_MCF7_Prog_5.txt"     1 MCF7_E2+Progesterone_3hr_Rep5
# [17] "GSM1669176_MCF7_R5020_6.txt"     MCF7_E2+R5020_3hr_Rep6
# [18] "GSM1669177_MCF7_Control_4.txt"  1 MCF7_E2_Rep4
# [19] "GSM1669178_MCF7_R5020_2.txt"     MCF7_E2+R5020_3hr_Rep2
# [20] "GSM1669179_MCF7_R5020_8.txt"     MCF7_E2+R5020_3hr_Rep8
# [21] "GSM1669180_MCF7_Prog_8.txt"     1 MCF7_E2+Progesterone_3hr_Rep8
# [22] "GSM1669181_MCF7_Prog_7.txt"     1 MCF7_E2+Progesterone_3hr_Rep7
# [23] "GSM1669182_MCF7_Prog_6.txt"     1 MCF7_E2+Progesterone_3hr_Rep6
# [24] "GSM1669183_MCF7_Control_2.txt"  1 MCF7_E2_Rep2
###################### THIS IS 8 REPLICATES OF THE SAME EXPERIMENT.
#exp1: control: [11]MCF7_E2_Rep1 other: [4]MCF7_E2+Progesterone_3hr_Rep1, [8]MCF7_E2+R5020_3hr_Rep1
#exp2: control: [24]MCF7_E2_Rep2 other: [15]MCF7_E2+Progesterone_3hr_Rep2, [19]MCF7_E2+R5020_3hr_Rep2
#exp3: control: [1]MCF7_E2_Rep3  other: [3]MCF7_E2+Progesterone_3hr_Rep3, [12]MCF7_E2+R5020_3hr_Rep3
#exp4: control: [18]MCF7_E2_Rep4 other: [6]MCF7_E2+Progesterone_3hr_Rep4, [7]MCF7_E2+R5020_3hr_Rep4
#exp5: control: [5]MCF7_E2_Rep5  other: [16]MCF7_E2+Progesterone_3hr_Rep5,[2]MCF7_E2+R5020_3hr_Rep5
#exp6: control: [9]MCF7_E2_Rep6  other: [23]MCF7_E2+Progesterone_3hr_Rep6, [17]MCF7_E2+R5020_3hr_Rep6
#exp7: control: [14]MCF7_E2_Rep7 other:[22]MCF7_E2+Progesterone_3hr_Rep7,[13]MCF7_E2+R5020_3hr_Rep7
#exp8: control: [10]MCF7_E2_Rep8 other: [21]MCF7_E2+Progesterone_3hr_Rep8, [20] MCF7_E2+R5020_3hr_Rep8

aa <- c("E2_1", "E2+R5020_1", "E2+Proges_1",
        "E2_2", "E2+R5020_2", "E2+Proges_2",
        "E2_3", "E2+R5020_3", "E2+Proges_3",
        "E2_4", "E2+R5020_4", "E2+Proges_4",
        "E2_5", "E2+R5020_5", "E2+Proges_5",
        "E2_6", "E2+R5020_6", "E2+Proges_6",
        "E2_7", "E2+R5020_7", "E2+Proges_7",
        "E2_8", "E2+R5020_8", "E2+Proges_8")
####################################################################################################
####################################################################################################
####################################################################################################

colnames(GSE76507.RNAseq)#6 G9 and PHF20 kd   ctrl: 3,5,6
#1"shG9a_NOE2"       
#2"shG9a_E2"          
#3"shNT_NOE2"         
#4"shNT_E2"          
#5"shNT_NoE2.rep1"  
#6"shNT_NOE2.rep2"    
#7"shNT_E2.rep1"     
#8"shNT_E2.rep2"     
#9"shPHF20_NOE2.rep1" 
#10"shPHF20_NOE2.rep2"
#11"shPHF20_E2.rep1"   
#12"shPHF20_E2.rep2"  
##########
#exp1: control: [1]shG9a_NOE2 other: [2]shG9a_E2
#exp2: control: [3]shNT_NOE2, other: [4]shNT_E2
#exp3: control: [5]shNT_NoE2.rep1, other: [7]shNT_E2.rep1
#exp4: control: [6]shNT_NoE2.rep2, other: [8]shNT_E2.rep2
#exp5: control: [9]shPHF20_NOE2.rep1, other: [11]shPHF20_E2.rep1
#exp6: control: [10]shPHF20_NOE2.rep2, other: [12]shPHF20_E2.rep2

aa <- c("shG9a_Veh", "shG9a_E2",
        "shNT_Veh_1", "shNT_E2_1",
        "shNT_Veh_2", "shNT_E2_2",
        "shNT_Veh_3", "shNT_E2_3",
        "shPHF20_Veh_1", "shPHF20_E2_1",
        "shPHF20_Veh_2", "shPHF20_E2_2")
####################################################################################################
####################################################################################################
####################################################################################################
colnames(GSE108308.RNAseq)      #R4K1 stapled peptide
#1 Veh_Rep1
#2 Veh_Rep2
#3 E2_Rep1
#4 E2_Rep2
#5 TOT_Rep1
#6 TOT_Rep2
#7 ET_Rep1
#8 ET_Rep2
#9 Peptide_Rep1
#10 Peptide_Rep2
#11 E_Peptide_Rep1
#12 E_Peptide_Rep2
#################
#exp1: control: [1]Veh_Rep1, other: [3]E2_Rep1 ,[5]TOT_Rep1,[7]ET_Rep1, [9]Peptide_Rep1, [11]E_Peptide_Rep1
#exp2: control: [3]Veh_Rep2, other: [4]E2_Rep2 ,[6]TOT_Rep2,[8]ET_Rep2,[10]Peptide_Rep2, [12]E_Peptide_Rep2

c("Veh_1", "E2_1", "Tamxfn_1", "E2+Tamxfn_1", "R4K1_1", "E2_R4K1_1",
  "Veh_2", "E2_2", "Tamxfn_2", "E2+Tamxfn_2", "R4K1_2", "E2_R4K1_1")
####################################################################################################
####################################################################################################
####################################################################################################
Exp.kd.Raw.Dataset.list.Common #is the list containing all the 7 datasets each as one entry in the same order as above. Only genes common between all datasets
#GSE73663.Inh.RNAseq  #1 Zeynep PaPe
#GSE86316.RNAseq      #2 MEN1 KO ##DEOPPED
#GSE55922.RNAseq      #3 BRD4 kd ##DROPPED
#GSE67295.RNAseq      #4 estradiol- and pro-inflammatory cytokine-treated
#GSE68358.RNAseq      #5 E2 or E2 plus something else
#GSE76507.RNAseq      #6 G9 and PHF kd
#GSE108308.RNAseq     #7 R4K1 stapled peptide
#Exp.kd.Batch.Dataset.matrix.Common.vgt1 #is the matrix of all experiments afer batch removal

#Create experiment list where each entry is a experiment matrix: for each expriment make a
#matrix that first column is the control experiment and the rest are others
DataSet.Experiment.List.Raw.vgt1   <- list()
DataSet.Experiment.List.Batch.vgt1 <- list()

#GSE73663RNA
#exp1 : control: 5"AZDVeh4h1" other : 1 "AZD334h1", 3 AZDE24h1
#exp2 : control: 6"AZDVeh4h2" other : 2"AZD334h2",4 AZDE24h2
#exp3 : control: 17 PP242Veh4h1 other: 13 PP242334h1 , 15 PP242E24h1
#exp4 : control: 18 PP242Veh4h2 other: 14 PP242334h2 ,16 PP242E24h2
#exp5 : control: 11 CtrlVeh4h1 other:  7 Ctrl334h1 , 9 CtrlE24h1
#exp6 : control: 12 CtrlVeh4h2 other:  8 Ctrl334h2 , 10 CtrlE24h2

DataSet.Experiment.List.Raw.vgt1[[1]] <- list()
DataSet.Experiment.List.Raw.vgt1[[1]][[1]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[1]][,c(5,1,3)]
DataSet.Experiment.List.Raw.vgt1[[1]][[2]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[1]][,c(6,2,4)]
DataSet.Experiment.List.Raw.vgt1[[1]][[3]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[1]][,c(17,13,15)]
DataSet.Experiment.List.Raw.vgt1[[1]][[4]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[1]][,c(18,14,16)]
DataSet.Experiment.List.Raw.vgt1[[1]][[5]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[1]][,c(11,7,9)]
DataSet.Experiment.List.Raw.vgt1[[1]][[6]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[1]][,c(12,8,10)]

colnames(DataSet.Experiment.List.Raw.vgt1[[1]][[1]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[1]][,c(5,1,3)])
colnames(DataSet.Experiment.List.Raw.vgt1[[1]][[2]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[1]][,c(6,2,4)])
colnames(DataSet.Experiment.List.Raw.vgt1[[1]][[3]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[1]][,c(17,13,15)])
colnames(DataSet.Experiment.List.Raw.vgt1[[1]][[4]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[1]][,c(18,14,16)])
colnames(DataSet.Experiment.List.Raw.vgt1[[1]][[5]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[1]][,c(11,7,9)])
colnames(DataSet.Experiment.List.Raw.vgt1[[1]][[6]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[1]][,c(12,8,10)])
#the version with batch removed
DataSet.Experiment.List.Batch.vgt1[[1]] <- list()
DataSet.Experiment.List.Batch.vgt1[[1]][[1]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[1]][,c(5,1,3)]
DataSet.Experiment.List.Batch.vgt1[[1]][[2]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[1]][,c(6,2,4)]
DataSet.Experiment.List.Batch.vgt1[[1]][[3]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[1]][,c(17,13,15)]
DataSet.Experiment.List.Batch.vgt1[[1]][[4]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[1]][,c(18,14,16)]
DataSet.Experiment.List.Batch.vgt1[[1]][[5]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[1]][,c(11,7,9)]
DataSet.Experiment.List.Batch.vgt1[[1]][[6]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[1]][,c(12,8,10)]
############
#mutualGSE86316RNA #2 Men1 kd ################FOR NOW DROPPING THIS
#exp1 : control: 1 Veh_ShCtrl_1        other: 3 E2_ShCtrl_1 
#exp2 : control: 2 Veh_ShCtrl_2        other: 4 E2_ShCtrl_2
#exp3 : control: 5 Veh_ShMEN1_1        other: 7 E2_ShMEN1_1
#exp4 : control: 6 Veh_ShMEN1_2        other: 8 E2_ShMEN1_2
# DataSet.Experiment.List.Raw.vgt1[[2]] <- list()
# DataSet.Experiment.List.Raw.vgt1[[2]][[1]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(1,3)]
# DataSet.Experiment.List.Raw.vgt1[[2]][[2]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(2,4)]
# DataSet.Experiment.List.Raw.vgt1[[2]][[3]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(5,7)]
# DataSet.Experiment.List.Raw.vgt1[[2]][[4]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(6,8)]
# 
# colnames(DataSet.Experiment.List.Raw.vgt1[[2]][[1]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(1,3)])
# colnames(DataSet.Experiment.List.Raw.vgt1[[2]][[2]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(2,4)])
# colnames(DataSet.Experiment.List.Raw.vgt1[[2]][[3]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(5,7)])
# colnames(DataSet.Experiment.List.Raw.vgt1[[2]][[4]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(6,8)])
# 
# #the version with batch removed
# DataSet.Experiment.List.Batch.vgt1[[2]] <- list()
# DataSet.Experiment.List.Batch.vgt1[[2]][[1]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[2]][,c(1,3)]
# DataSet.Experiment.List.Batch.vgt1[[2]][[2]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[2]][,c(2,4)]
# DataSet.Experiment.List.Batch.vgt1[[2]][[3]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[2]][,c(5,7)]
# DataSet.Experiment.List.Batch.vgt1[[2]][[4]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[2]][,c(6,8)]
###############
#mutualGSE53532RNA) #3 CTCF kd     #THis is ommited from dataset
#"siNT_E0h"   
#"siNT_E3h"   
#"siCTCF_E0h" 
#"siCTCF_E3h"
#############
#exp1: control: 1 siNT_E0h      other:2 siNT_E3h
#exp2: control: 3 siCTCF_E0h    other:4 siCTCF_E3h
#DataSetExpRaw[[3]] <- list()
#DataSetExpRaw[[3]][[1]] <- mutualGSE53532RNA[,c(1,2)]
#DataSetExpRaw[[3]][[2]] <- mutualGSE53532RNA[,c(3,4)]

#colnames(DataSetExpRaw[[3]][[1]]) <- colnames(mutualGSE53532RNA[,c(1,2)])
#colnames(DataSetExpRaw[[3]][[2]]) <- colnames(mutualGSE53532RNA[,c(3,4)])

#the version with batch removed
#DataSetExpBatch[[3]] <- list()
#DataSetExpBatch[[3]][[1]] <- mutualAll7KDdatasetsBatch[,(c(1,2)+26)]
#DataSetExpBatch[[3]][[2]] <- mutualAll7KDdatasetsBatch[,(c(3,4)+26)]
#####################
#mutualGSE55922RNA #3 BRD4 kd     #FOR NOW DROPPING THIS
#exp1: control: 1 brd4.e2a other: 3 brd4pluse2a
#exp2: control: 2 brd4.e2b other: 4 brd4pluse2b
#exp3: control: 5 cont.e2a other: 7 contpluse2a
#exp4: control: 6 cont.e2b other: 8 contpluse2b
#exp5: control: 9 jq1.e2a other: 11 jq1pluse2a
#exp6: control: 10 jq1.e2b other: 12 jq1pluse2b
# DataSet.Experiment.List.Raw.vgt1[[3]] <- list()
# DataSet.Experiment.List.Raw.vgt1[[3]][[1]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(1,3)]
# DataSet.Experiment.List.Raw.vgt1[[3]][[2]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(2,4)]
# DataSet.Experiment.List.Raw.vgt1[[3]][[3]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(5,7)]
# DataSet.Experiment.List.Raw.vgt1[[3]][[4]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(6,8)]
# DataSet.Experiment.List.Raw.vgt1[[3]][[5]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(9,11)]
# DataSet.Experiment.List.Raw.vgt1[[3]][[6]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(10,12)]
# 
# colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[1]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(1,3)])
# colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[2]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(2,4)])
# colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[3]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(5,7)])
# colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[4]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(6,8)])
# colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[5]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(9,11)])
# colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[6]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(10,12)])
# 
# #the version with batch removed
# DataSet.Experiment.List.Batch.vgt1[[3]] <- list()
# DataSet.Experiment.List.Batch.vgt1[[3]][[1]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(1,3)]
# DataSet.Experiment.List.Batch.vgt1[[3]][[2]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(2,4)]
# DataSet.Experiment.List.Batch.vgt1[[3]][[3]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(5,7)]
# DataSet.Experiment.List.Batch.vgt1[[3]][[4]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(6,8)]
# DataSet.Experiment.List.Batch.vgt1[[3]][[5]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(9,11)]
# DataSet.Experiment.List.Batch.vgt1[[3]][[6]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(10,12)]
########42
#mutualGSE67295RNA #5 #E2 and other treatments 
############################
#exp1: control: 1 "Sample19.RNASeq.Veh.Expt1.Rep1.fastq.gz" other 3 Sample21.RNASeq.E2.Rep1.fastq.gz, 5 Sample23.RNASeq.E2.ICI.Rep1.fastq.gz, 7 Sample25.RNASeq.IL1b.Rep1.fastq.gz, 9 Sample27.RNASeq.IL1b.ICI.Rep1.fastq.gz
#exp2: control: 2 Sample20.RNASeq.Veh.Expt1.Rep2.fastq.gz   other 4 Sample22.RNASeq.E2.Rep2.fastq.gz, 6 Sample24.RNASeq.E2.ICI.Rep2.fastq.gz, 8 Sample26.RNASeq.IL1b.Rep2.fastq.gz, 10 Sample28.RNASeq.IL1b.ICI.Rep2.fastq.gz
#exp3: control: 11 Sample29.RNASeq.Veh.Rep1.fastq.gz         other 13 Sample31.RNASeq.TNF.Rep1.fastq.gz, 15 Sample33.RNASeq.TNF.ICI.Rep1.fastq.gz
#exp4: control: 12 Sample30.RNASeq.Veh.Rep2.fastq.gz         other 14 Sample32.RNASeq.TNF.Rep2.fastq.gz, 16 Sample34.RNASeq.TNF.ICI.Rep2.fastq.gz
#exp5: control: 17 Sample35.RNASeq.Veh.Expt3.Rep1.fastq.gz   other 19 Sample37.RNASeq.E2.Expt3.Rep1.fastq.gz, 21 Sample39.RNASeqE2.TOT.Expt3.Rep1.fastq.gz,23 Sample41.RNASeq.E2.TOT.IL1b.Expt3.Rep1.fastq.gz, 25 Sample43.RNASeq.E2.TOT.TNF.Expt3.Rep1.fastq.gz
#exp6: control: 18 Sample36.RNASeq.Veh.Expt3.Rep2.fastq.gz   other 20 Sample38.RNASeq.E2.Expt3.Rep2.fastq.gz, 22 Sample40.RNASeq.E2.TOT.Expt3.Rep2.fastq.gz,24 Sample42.RNASeq.E2.TOT.IL1b.Expt3.Rep2.fastq.gz, 26 Sample44.RNASeq.E2.TOT.TNF.Expt3.Rep2.fastq.gz

DataSet.Experiment.List.Raw.vgt1[[2]] <- list()
DataSet.Experiment.List.Raw.vgt1[[2]][[1]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(1,3,5,7,9)]
DataSet.Experiment.List.Raw.vgt1[[2]][[2]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(2,4,6,8,10)]
DataSet.Experiment.List.Raw.vgt1[[2]][[3]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(11,13,15)]
DataSet.Experiment.List.Raw.vgt1[[2]][[4]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(12,14,16)]
DataSet.Experiment.List.Raw.vgt1[[2]][[5]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(17,19,21,23,25)]
DataSet.Experiment.List.Raw.vgt1[[2]][[6]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(18,20,22,24,26)]

colnames(DataSet.Experiment.List.Raw.vgt1[[2]][[1]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(1,3,5,7,9)])
colnames(DataSet.Experiment.List.Raw.vgt1[[2]][[2]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(2,4,6,8,10)])
colnames(DataSet.Experiment.List.Raw.vgt1[[2]][[3]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(11,13,15)])
colnames(DataSet.Experiment.List.Raw.vgt1[[2]][[4]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(12,14,16)])
colnames(DataSet.Experiment.List.Raw.vgt1[[2]][[5]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(17,19,21,23,25)])
colnames(DataSet.Experiment.List.Raw.vgt1[[2]][[6]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[2]][,c(18,20,22,24,26)])


#the version with batch removed
DataSet.Experiment.List.Batch.vgt1[[2]] <- list()
DataSet.Experiment.List.Batch.vgt1[[2]][[1]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[2]][,c(1,3,5,7,9)]
DataSet.Experiment.List.Batch.vgt1[[2]][[2]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[2]][,c(2,4,6,8,10)]
DataSet.Experiment.List.Batch.vgt1[[2]][[3]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[2]][,c(11,13,15)]
DataSet.Experiment.List.Batch.vgt1[[2]][[4]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[2]][,c(12,14,16)]
DataSet.Experiment.List.Batch.vgt1[[2]][[5]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[2]][,c(17,19,21,23,25)]
DataSet.Experiment.List.Batch.vgt1[[2]][[6]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[2]][,c(18,20,22,24,26)]
##################
#mutualGSE68358RNA#6 #E2 or E2 plus something 
######################
#exp1: control: [11]MCF7_E2_Rep1 other: [4]MCF7_E2+Progesterone_3hr_Rep1, [8]MCF7_E2+R5020_3hr_Rep1
#exp2: control: [24]MCF7_E2_Rep2 other: [15]MCF7_E2+Progesterone_3hr_Rep2, [19]MCF7_E2+R5020_3hr_Rep2
#exp3: control: [1]MCF7_E2_Rep3  other: [3]MCF7_E2+Progesterone_3hr_Rep3, [12]MCF7_E2+R5020_3hr_Rep3
#exp4: control: [18]MCF7_E2_Rep4 other: [6]MCF7_E2+Progesterone_3hr_Rep4, [7]MCF7_E2+R5020_3hr_Rep4
#exp5: control: [5]MCF7_E2_Rep5  other: [16]MCF7_E2+Progesterone_3hr_Rep5,[2]MCF7_E2+R5020_3hr_Rep5
#exp6: control: [9]MCF7_E2_Rep6  other: [23]MCF7_E2+Progesterone_3hr_Rep6, [17]MCF7_E2+R5020_3hr_Rep6
#exp7: control: [14]MCF7_E2_Rep7 other:[22]MCF7_E2+Progesterone_3hr_Rep7,[13]MCF7_E2+R5020_3hr_Rep7
#exp8: control: [10]MCF7_E2_Rep8 other: [21]MCF7_E2+Progesterone_3hr_Rep8, [20] MCF7_E2+R5020_3hr_Rep8
DataSet.Experiment.List.Raw.vgt1[[3]] <- list()
DataSet.Experiment.List.Raw.vgt1[[3]][[1]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(11,4,8)]
DataSet.Experiment.List.Raw.vgt1[[3]][[2]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(24,15,19)]
DataSet.Experiment.List.Raw.vgt1[[3]][[3]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(1,3,12)]
DataSet.Experiment.List.Raw.vgt1[[3]][[4]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(18,6,7)]
DataSet.Experiment.List.Raw.vgt1[[3]][[5]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(5,16,2)]
DataSet.Experiment.List.Raw.vgt1[[3]][[6]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(9,23,17)]
DataSet.Experiment.List.Raw.vgt1[[3]][[7]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(14,22,13)]
DataSet.Experiment.List.Raw.vgt1[[3]][[8]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(10,21,20)]

colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[1]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(11,4,8)])
colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[2]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(24,15,19)])
colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[3]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(1,3,12)])
colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[4]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(18,6,7)])
colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[5]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(5,16,2)])
colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[6]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(9,23,17)])
colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[7]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(14,22,13)])
colnames(DataSet.Experiment.List.Raw.vgt1[[3]][[8]]) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[3]][,c(10,21,20)])

#the version with batch removed
DataSet.Experiment.List.Batch.vgt1[[3]] <- list()
DataSet.Experiment.List.Batch.vgt1[[3]][[1]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(11,4,8)]
DataSet.Experiment.List.Batch.vgt1[[3]][[2]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(24,15,19)]
DataSet.Experiment.List.Batch.vgt1[[3]][[3]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(1,3,12)]
DataSet.Experiment.List.Batch.vgt1[[3]][[4]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(18,6,7)]
DataSet.Experiment.List.Batch.vgt1[[3]][[5]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(5,16,2)]
DataSet.Experiment.List.Batch.vgt1[[3]][[6]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(9,23,17)]
DataSet.Experiment.List.Batch.vgt1[[3]][[7]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(14,22,13)]
DataSet.Experiment.List.Batch.vgt1[[3]][[8]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[3]][,c(10,21,20)]
###########################92
#mutualGSE76507RNA)#7 G9 and PHF20 kd
##########
#exp1: control: [1]shG9a_NOE2 other: [2]shG9a_E2
#exp2: control: [3]shNT_NOE2, other: [4]shNT_E2
#exp3: control: [5]shNT_NoE2.rep1, other: [7]shNT_E2.rep1
#exp4: control: [6]shNT_NoE2.rep2, other: [8]shNT_E2.rep2
#exp5: control: [9]shPHF20_NOE2.rep1, other: [11]shPHF20_E2.rep1
#exp6: control: [10]shPHF20_NOE2.rep2, other: [12]shPHF20_E2.rep2
DataSet.Experiment.List.Raw.vgt1[[4]] <- list()
DataSet.Experiment.List.Raw.vgt1[[4]][[1]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[4]][,c(1,2)]
DataSet.Experiment.List.Raw.vgt1[[4]][[2]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[4]][,c(3,4)]
DataSet.Experiment.List.Raw.vgt1[[4]][[3]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[4]][,c(5,7)]
DataSet.Experiment.List.Raw.vgt1[[4]][[4]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[4]][,c(6,8)]
DataSet.Experiment.List.Raw.vgt1[[4]][[5]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[4]][,c(9,11)]
DataSet.Experiment.List.Raw.vgt1[[4]][[6]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[4]][,c(10,12)]

colnames(DataSet.Experiment.List.Raw.vgt1[[4]][[1]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[4]][,c(1,2)])
colnames(DataSet.Experiment.List.Raw.vgt1[[4]][[2]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[4]][,c(3,4)])
colnames(DataSet.Experiment.List.Raw.vgt1[[4]][[3]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[4]][,c(5,7)])
colnames(DataSet.Experiment.List.Raw.vgt1[[4]][[4]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[4]][,c(6,8)])
colnames(DataSet.Experiment.List.Raw.vgt1[[4]][[5]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[4]][,c(9,11)])
colnames(DataSet.Experiment.List.Raw.vgt1[[4]][[6]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[4]][,c(10,12)])

#the version with batch removed
DataSet.Experiment.List.Batch.vgt1[[4]] <- list()
DataSet.Experiment.List.Batch.vgt1[[4]][[1]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[4]][,c(1,2)]
DataSet.Experiment.List.Batch.vgt1[[4]][[2]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[4]][,c(3,4)]
DataSet.Experiment.List.Batch.vgt1[[4]][[3]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[4]][,c(5,7)]
DataSet.Experiment.List.Batch.vgt1[[4]][[4]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[4]][,c(6,8)]
DataSet.Experiment.List.Batch.vgt1[[4]][[5]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[4]][,c(9,11)]
DataSet.Experiment.List.Batch.vgt1[[4]][[6]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[4]][,c(10,12)]
###########################
#colnames(GSE108308.RNAseq)     7 #R4K1 stapled peptide
#1 Veh_Rep1
#2 Veh_Rep2
#3 E2_Rep1
#4 E2_Rep2
#5 TOT_Rep1
#6 TOT_Rep2
#7 ET_Rep1
#8 ET_Rep2
#9 Peptide_Rep1
#10 Peptide_Rep2
#11 E_Peptide_Rep1
#12 E_Peptide_Rep2
#################
#exp1: control: [1]Veh_Rep1, other: [3]E2_Rep1 ,[5]TOT_Rep1,[7]ET_Rep1, [9]Peptide_Rep1, [11]E_Peptide_Rep1
#exp2: control: [3]Veh_Rep2, other: [4]E2_Rep2 ,[6]TOT_Rep2,[8]ET_Rep2,[10]Peptide_Rep2, [12]E_Peptide_Rep2
DataSet.Experiment.List.Raw.vgt1[[5]] <- list()
DataSet.Experiment.List.Raw.vgt1[[5]][[1]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[5]][,c(1,3,5,7,9,11)]
DataSet.Experiment.List.Raw.vgt1[[5]][[2]] <- Exp.kd.Raw.Dataset.list.Common.vgt1[[5]][,c(2,4,6,8,10,12)]
#
colnames(DataSet.Experiment.List.Raw.vgt1[[5]][[1]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[5]][,c(1,3,5,7,9,11)])
colnames(DataSet.Experiment.List.Raw.vgt1[[5]][[2]] ) <- colnames(Exp.kd.Raw.Dataset.list.Common.vgt1[[5]][,c(2,4,6,8,10,12)])
#
DataSet.Experiment.List.Batch.vgt1[[5]] <- list()
DataSet.Experiment.List.Batch.vgt1[[5]][[1]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[5]][,c(1,3,5,7,9,11)]
DataSet.Experiment.List.Batch.vgt1[[5]][[2]] <- Exp.kd.Batch.Dataset.list.Common.vgt1[[5]][,c(2,4,6,8,10,12)]

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#put all experiments next to each other and note the number of each dataset, name of each experiment, and control experiment
a <- list()
for(i in 1:length(DataSet.Experiment.List.Raw.vgt1)){
  a[[i]] <- do.call(cbind, DataSet.Experiment.List.Raw.vgt1[[i]])
}
DataSet.Experiment.matrix.Raw.vgt1 <- do.call(cbind, a)
###
a <- list()
for(i in 1:length(DataSet.Experiment.List.Batch.vgt1)){
  a[[i]] <- do.call(cbind, DataSet.Experiment.List.Batch.vgt1[[i]])
}
DataSet.Experiment.matrix.Batch.vgt1 <- do.call(cbind, a)

#######################################################################################################################################
#create metadata: 112 rows, 2 columns, first column "DatasetNo", second column "expNo" where #expNo == 1 is the control exp

DataSet.Experiment.Metadata <- data.frame(DatasetNo = numeric(ncol(DataSet.Experiment.matrix.Batch.vgt1)), expNo = numeric(ncol(DataSet.Experiment.matrix.Batch.vgt1)))
aa <- unlist(lapply(Exp.kd.Raw.Dataset.list.Common, ncol))
a <- numeric(0)
for(i in 1:length(aa)){
  a <- c(a, rep(i,aa[i]))
}
DataSet.Experiment.Metadata$DatasetNo = a
#

a <- numeric(0)
for(i in 1:length(DataSet.Experiment.List.Raw.vgt1)){
  aa <- unlist(lapply(DataSet.Experiment.List.Raw.vgt1[[i]],ncol))
  for(j in 1:length(aa)){
    a = c(a, c(1:aa[j]))
  }
}

DataSet.Experiment.Metadata$expNo = a
rownames(DataSet.Experiment.Metadata) <- colnames(DataSet.Experiment.matrix.Raw.vgt1)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#Construct fold change matrix, format is the same as 
#DataSet.Experiment.matrix.Raw.vgt1, but 
#each entry is log fold change from the control experiment
DataSet.Experiment.List.Raw.vgt1.logfold <- DataSet.Experiment.List.Raw.vgt1
DataSet.Experiment.List.Batch.vgt1.logfold <- DataSet.Experiment.List.Batch.vgt1
for(i in 1:length(DataSet.Experiment.List.Raw.vgt1.logfold)){
  for(j in 1:length(DataSet.Experiment.List.Raw.vgt1.logfold[[i]])){
    #adding one to all values
    DataSet.Experiment.List.Raw.vgt1.logfold[[i]][[j]] <- DataSet.Experiment.List.Raw.vgt1.logfold[[i]][[j]] + 1
    #dividing by the (control experiment)
    DataSet.Experiment.List.Raw.vgt1.logfold[[i]][[j]] <- DataSet.Experiment.List.Raw.vgt1.logfold[[i]][[j]]/DataSet.Experiment.List.Raw.vgt1.logfold[[i]][[j]][,1]
    #computing log2(x)
    DataSet.Experiment.List.Raw.vgt1.logfold[[i]][[j]] <- log2(DataSet.Experiment.List.Raw.vgt1.logfold[[i]][[j]])
    #setting control column equal to one
    #DataSet.Experiment.List.Raw.vgt1.logfold[[i]][[j]][,1] <- matrix(1L, nrow = nrow(DataSet.Experiment.List.Raw.vgt1.logfold[[i]][[j]]), ncol = 1)
    #setting negative values to zero
    DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]][DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]] < 0] <- 0
    #adding one to all values
    DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]] <- DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]] + 1
    
    #dividing by the (control experiment)
    DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]] <- DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]]/DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]][,1]
    #computing log2(x)
   # print("beforelog")
   # print(sum(DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]][,1] == 1))
    DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]] <- log2(DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]])
   # print("afterlog")
    #print(sum(DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]][,1] == 0))
    #setting control column equal to one
    #DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]][,1] <- matrix(1L, nrow = nrow(DataSet.Experiment.List.Batch.vgt1.logfold[[i]][[j]]), ncol = 1)
  }
}

a <- list()
for(i in 1:length(DataSet.Experiment.List.Raw.vgt1.logfold)){
  a[[i]] <- do.call(cbind, DataSet.Experiment.List.Raw.vgt1.logfold[[i]])
}
DataSet.Experiment.Matrix.Raw.vgt1.logfold <- do.call(cbind, a)

a <- list()
for(i in 1:length(DataSet.Experiment.List.Batch.vgt1.logfold)){
  a[[i]] <- do.call(cbind, DataSet.Experiment.List.Batch.vgt1.logfold[[i]])
}
DataSet.Experiment.Matrix.Batch.vgt1.logfold <- do.call(cbind, a)

boxplot.matrix(DataSet.Experiment.Matrix.Batch.vgt1.logfold,main="LFC batch removed all genes", xlab = "",xaxt = "n" )
boxplot.matrix(DataSet.Experiment.Matrix.Raw.vgt1.logfold,main="LFC raw removed all genes", xlab = "", xaxt = "n" )

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

#######################################################################################################################################
###########################################   FULL   ############################################################################
#######################################################################################################################################
#FULL DATASETS:
DataSet.Experiment.matrix.Raw.vgt1
DataSet.Experiment.Matrix.Raw.vgt1.logfold
DataSet.Experiment.matrix.Batch.vgt1
DataSet.Experiment.Matrix.Batch.vgt1.logfold
DataSet.Experiment.matrix.Names <- c(c("MAPK_Veh_1", "MAPK_PaPE_1", "MAPK_E2_1",
                                       "MAPK_Veh_2", "MAPK_PaPE_2", "MAPK_E2_2",
                                       "mTOR_Veh_1", "mTOR_PaPE_1", "mTOR_E2_1",
                                       "mTOR_Veh_2", "mTOR_PaPE_2", "mTOR_E2_2",
                                       "Veh_Veh_1", "Veh_PaPE_1", "Veh_E2_1",
                                       "Veh_Veh_2", "Veh_PaPE_2", "Veh_E2_2"), 
                                     c("Veh_1_1", "E2_1_1", "Fulv+E2_1_1",  "Il1b_1_1", "Fulv+Il1b_1_1",
                                       "Veh_1_2", "E2_1_2", "Fulv+E2_1_2",  "Il1b_1_2", "Fulv+Il1b_1_2",
                                       "Veh_2_1", "TNFa_2_1", "Fulv+TNFa_2_1",
                                       "Veh_2_2", "TNFa_2_2", "Fulv+TNFa_2_2",
                                       "Veh_3_1",  "E2_3_1", "Tamxfn+E2_3_1", "Tamxfn+E2+IL1b_3_1", "Tamxfn+E2+TNFa_3_1",
                                       "Veh_3_2",  "E2_3_2", "Tamxfn+E2_3_2", "Tamxfn+E2+IL1b_3_2", "Tamxfn+E2+TNFa_3_2"),
                                     c("E2_1", "E2+R5020_1", "E2+Proges_1",
                                       "E2_2", "E2+R5020_2", "E2+Proges_2",
                                       "E2_3", "E2+R5020_3", "E2+Proges_3",
                                       "E2_4", "E2+R5020_4", "E2+Proges_4",
                                       "E2_5", "E2+R5020_5", "E2+Proges_5",
                                       "E2_6", "E2+R5020_6", "E2+Proges_6",
                                       "E2_7", "E2+R5020_7", "E2+Proges_7",
                                       "E2_8", "E2+R5020_8", "E2+Proges_8"),
                                     c("shG9a_Veh", "shG9a_E2",
                                       "shNT_Veh_1", "shNT_E2_1",
                                       "shNT_Veh_2", "shNT_E2_2",
                                       "shNT_Veh_3", "shNT_E2_3",
                                       "shPHF20_Veh_1", "shPHF20_E2_1",
                                       "shPHF20_Veh_2", "shPHF20_E2_2"),
                                     c("Veh_1", "E2_1", "Tamxfn_1", "E2+Tamxfn_1", "R4K1_1", "E2_R4K1_1",
                                       "Veh_2", "E2_2", "Tamxfn_2", "E2+Tamxfn_2", "R4K1_2", "E2_R4K1_1"))
#######################################################################################################################################
##########################################   ER ONLY   ################################################################################
#######################################################################################################################################
#datasets for genes associated with an ER chip peak only:
aa = which(rownames(DataSet.Experiment.matrix.Batch.vgt1) %in% Genes.Associated.REMAP.ER.Entrez)
DataSet.Experiment.matrix.Raw.vgt1.ER <- DataSet.Experiment.matrix.Raw.vgt1[aa,]
DataSet.Experiment.matrix.Raw.vgt1.LFC.ER <- DataSet.Experiment.Matrix.Raw.vgt1.logfold[aa,]
DataSet.Experiment.matrix.Batch.vgt1.ER <- DataSet.Experiment.matrix.Batch.vgt1[aa,]
DataSet.Experiment.matrix.Batch.vgt1.LFC.ER <- DataSet.Experiment.Matrix.Batch.vgt1.logfold[aa,]
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
boxplot.matrix(DataSet.Experiment.matrix.Batch.vgt1.LFC.ER,main="LFC batch removed ER genes", xlab = "", las = 2 )
boxplot.matrix(DataSet.Experiment.matrix.Raw.vgt1.LFC.ER,main="LFC raw removed ER genes", xlab = "", las = 2 )

