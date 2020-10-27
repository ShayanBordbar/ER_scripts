# e_RNA transcription prediction
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")

#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         Libraries         ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
library(ChIPseeker)
library(GenomicRanges)
source("https://bioconductor.org/biocLite.R")
biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("BSgenome.Hsapiens.UCSC.hg38")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         Functions         ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
bedtools_intersect <- function(bedfile_names, wa=F, wb=F, loj=F, wo=F, wao=F,
                              u=F, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                              S=F, sorted=F, merge_after_distance=0, output_filename, read_output=T){
  # bedfile_names: a character vector containing address of bed files, has to be at least of 
  #  length 2. the first entry will be treated as A and the rest of entries as B
  # wa logical, if True: Write the original entry in A for each overlap.
  # wb logical, if True: Write the original entry in B for each overlap. Useful for knowing
  #  what A overlaps. Restricted by f and r.
  # loj logical, if True: Perform a “left outer join”. That is, for each feature in A report each
  #  overlap with B. If no overlaps are found, report a NULL feature for B.
  # wo logical, if True: Write the original A and B entries plus the number of base pairs of
  #  overlap between the two features. Only A features with overlap are reported. Restricted by -f and -r
  # wao logical, if True: Write the original A and B entries plus the number of base pairs of
  #  overlap between the two features. However, A features w/o overlap are also reported with a
  #  NULL B feature and overlap = 0. Restricted by -f and -r.
  # u logical, if True: Write original A entry once if any overlaps found in B. In other words,
  #  just report the fact at least one overlap was found in B. Restricted by -f and -r
  # c logical, if True: For each entry in A, report the number of hits in B while 
  #  restricting to -f. Reports 0 for A entries that have no overlap with B. Restricted by -f and -r.
  # v logical, if True:Only report those entries in A that have no overlap in B. Restricted by -f and -r.
  # f numeric :  Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).
  # .F numeric : Minimum overlap required as a fraction of B. Default is 1E-9 (i.e., 1bp).
  # r logical, Require that the fraction of overlap be reciprocal for A and B. In other words,
  #  if -f is 0.90 and -r is used, this requires that B overlap at least 90% of A and that A also
  #  overlaps at least 90% of B.
  # e logical: Require that the minimum fraction be satisfied for A _OR_ B. In other words, 
  #  if -e is used with -f 0.90 and -F 0.10 this requires that either 90% of A is covered OR
  #  10% of B is covered. Without -e, both fractions would have to be satisfied.
  # s logical: Force “strandedness”. That is, only report hits in B that overlap A on the same
  #  strand. By default, overlaps are reported without respect to strand.
  # S logical : Require different strandedness. That is, only report hits in B that overlap A on 
  #  the _opposite_ strand. By default, overlaps are reported without respect to strand.
  # sorted logical : For very large B files, invoke a “sweeping” algorithm that requires
  #  position-sorted (e.g., sort -k1,1 -k2,2n for BED files) input. When using -sorted, memory usage
  #  remains low even for very large files.
  # merge_after_distance : integer, if greater than 0, it merges the peaks that are in the specified
  #  distance and adds a column refering to their counts
  # output_filename : name or address of output file
  # read_output: logical, if True reads and returns the output as a GRanges object using 
  #  readPeakFile function of ChIPseeker library
  # By default, if an overlap is found, bedtools intersect reports the shared interval between
  #  the two overlapping features.
  
  
  stopifnot(is.character(bedfile_names),
            length(bedfile_names) >= 2,
            #(sum(as.integer(c(wa, wb, loj, wo, wao, u, v))) <= 1),
            (sum(as.integer(c(r, e))) <= 1), 
            (sum(as.integer(c(s, S))) <= 1), 
            (f > 0 & .F > 0) 
            )
  
  output_filename_sep <- unlist(strsplit(output_filename, split = "\\."))
  if(length(output_filename_sep) == 1){
    output_filename_pref <- output_filename
    output_filename_full <- paste0(output_filename, ".bed")
  }else if(length(output_filename_sep) > 1 & output_filename_sep[length(output_filename_sep)] == "bed"){
    output_filename_pref <- paste0(output_filename_sep[1:(length(output_filename_sep) - 1)], collapse = "")
    output_filename_full <- paste0(output_filename_pref, ".bed")
  }else{
    output_filename_pref <- paste0(output_filename_sep, collapse = "")
    output_filename_full <- paste0(output_filename_pref, ".bed")
  }
  
  for(cur_file in 1:length(bedfile_names)){
    print(paste0("number of peaks in file ", bedfile_names[cur_file], " :"))
    system(paste0("cat ", bedfile_names[cur_file], " | wc -l"))
  }
  sys_command <- paste0("bedtools intersect ",
                        rep("-wa ", as.integer(wa)),
                        rep("-wb ", as.integer(wb)),
                        rep("-loj ", as.integer(loj)),
                        rep("-wo ", as.integer(wo)),
                        rep("-wao ", as.integer(wao)),
                        rep("-u ", as.integer(u)),
                        rep("-c ", as.integer(c)),
                        rep("-v ", as.integer(v)), 
                        paste0("-f ", f, " "),
                        paste0("-F ", .F, " "),
                        rep("-r ", as.integer(r)),
                        rep("-e ", as.integer(e)), 
                        rep("-s ", as.integer(s)),
                        rep("-S ", as.integer(S)),
                        rep("-sorted ", as.integer(sorted)),
                        " -a ", bedfile_names[1],
                        " -b ", paste(bedfile_names[2:length(bedfile_names)],
                                      collapse = " "), " ",
                        " | sortBed -i stdin > ", output_filename_full
                        )
  # debug
  print(sys_command)
  # 
  system(sys_command)
  print("number of peaks without merging closeby peaks: ")
  system(paste0("cat ", output_filename, " | wc -l"))
  

  if(merge_after_distance > 0){
    merged_name <- paste0(output_filename_pref, "_merged_", merge_after_distance, "bp.bed")
    sys_command_merge <- paste0("bedtools merge -i ",
                                output_filename,
                                " -d ", merge_after_distance, 
                                " -c 1 -o count > ",
                                merged_name
                                )
    system(sys_command_merge)
    
    print(paste0("number of peaks after merging peaks closer than ", as.character(merge_after_distance), " bp: "))
    system(paste0("cat ", merged_name, " | wc -l"))
    system(paste0("mv ", merged_name, " ", output_filename_full))
    # system(paste0("rm ", merged_name))
    #output_filename_full <- merged_name
  }
  if(read_output){
    my_result <- readPeakFile(output_filename_full, as = "GRanges")
    return(my_result)
  }
}
#########################################################################################################
#########################################################################################################
# example
aa_fname <- list.files("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/Li2013/", pattern = "*filtered_5.bed")[1:3]
aa_fname <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/Li2013/", aa_fname)

aa <- bedtools_intersect(bedfile_names = aa_fname, wa=F, wb=F, loj=F, wo=F, wao=F,
                         u=F, c=F, v=T, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                         S=F, sorted=F, merge_after_distance=100, 
                         output_filename="test_output.bed", read_output = T)


#########################################################################################################
#########################################################################################################
bedtools_merger <- function(bedfile_names, 
                            merge_dist,
                            col_number=1,
                            o="count",
                            output_filename,
                            read_output){
  # is a character vector containing the address of bedfiles to be merged
  # merge_dist : is an integer indicating the distance that peaks closer than should merge
  # col_number: is the -c option to bedtools merge tool. The -c option allows one to specify a column or columns in the input that you wish to summarize. it should be an integer vector of length one or more
  # The -o option defines the operation(s) that you wish to apply to each column listed for the -c option. For example,
  #  to count the number of overlapping intervals that led to each of the new “merged” intervals, one will “count” 
  #  the first column (though the second, third, fourth, etc. would work just fine as well).
  # output_filename is the name or full address of the output file
  # read_output : logical, if True reads and returns the output as a GRanges object using 
  #  readPeakFile function of ChIPseeker library
  stopifnot(is.numeric(merge_dist),
            length(merge_dist) == 1, 
            is.character(bedfile_names), 
            length(bedfile_names) >= 1,
            is.numeric(col_number),
            is.character(o), 
            is.character(output_filename),
            length(output_filename) == 1)
  merge_dist <- as.integer(merge_dist)
  col_number <- as.integer(col_number)
  output_filename_sep <- unlist(strsplit(output_filename, split = "\\."))
  if(length(output_filename_sep) == 1){
    output_filename_pref <- output_filename
    output_filename_full <- paste0(output_filename, ".bed")
  }else if(length(output_filename_sep) > 1 & output_filename_sep[length(output_filename_sep)] == "bed"){
    output_filename_pref <- paste0(output_filename_sep[1:(length(output_filename_sep) - 1)], collapse = "")
    output_filename_full <- paste0(output_filename_pref, ".bed")
  }else{
    output_filename_pref <- paste0(output_filename_sep, collapse = "")
    output_filename_full <- paste0(output_filename_pref, ".bed")
  }
  temp_list <- list()
  for(cur_file in 1:length(bedfile_names)){
    print(paste0("number of peaks in file ", bedfile_names[cur_file], " :"))
    system(paste0("cat ", bedfile_names[cur_file], " | wc -l"))
    temp_list[[cur_file]] <- readPeakFile(bedfile_names[cur_file],as = "data.frame")
  }
  names(temp_list) <- bedfile_names
  field_nu <- unlist(lapply(temp_list, ncol))
  bedfile_names_nobed <- unlist(strsplit(bedfile_names, split = ".bed"))
  if(length(unique(field_nu)) > 1){
    min_fields <- min(field_nu)
    bad_beds <- bedfile_names[field_nu > min_fields]
    bad_beds_nobed <- bedfile_names_nobed[field_nu > min_fields]
    bad_beds_list <- temp_list[field_nu > min_fields]
    print(paste0("Input bed files ", paste0(bad_beds, collapse = " "), " have more than ",
                 min_fields, " fields. Another version of them containing only ",
                 min_fields, " fields is created and used."))
    options(scipen=999)
    for(cur_b in 1:length(bad_beds_list)){
      write.table(bad_beds_list[[cur_b]][, 1: min_fields], 
                  file=paste0(bad_beds_nobed[cur_b], "_",min_fields,"_fields", ".bed"),
                  quote=F, sep="\t", row.names=F, col.names=F)
    }
    options(scipen=0)
    bedfile_names[field_nu > min_fields] <- paste0(bad_beds_nobed, "_",min_fields,"_fields", ".bed")
  }
  # first concatenate and sort the bed files if there are multiple of them or if there is one and its not sorted
  sys_command_concat_sort <- paste0("cat ", 
                                    paste(bedfile_names, collapse = " "),
                                    " | sortBed -i stdin > canc_temp.bed")
  print(sys_command_concat_sort)
  system(sys_command_concat_sort)
  # now merge using the provided arguments
  sys_command_merger <- paste0("bedtools merge -i canc_temp.bed -d ",
                               merge_dist, 
                               rep(" -c ", as.integer(length(col_number) > 0)), col_number,
                               rep(" -o ", as.integer(length(o) > 0)), o, 
                               " > ", output_filename_full)
  print(sys_command_merger)
  system(sys_command_merger)
  system("rm canc_temp.bed")
  print(paste0("number of peaks after merging peaks closer than ", as.character(merge_dist), " bp: "))
  system(paste0("cat ", output_filename_full, " | wc -l"))
  
  if(read_output){
    my_result <- readPeakFile(output_filename_full, as = "GRanges")
    return(my_result)
  }
}
#########################################################################################################
aa_fname <- list.files("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/Li2013/", pattern = "*filtered_5.bed")[1:3]
aa_fname <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/Li2013/", aa_fname)

system(paste0("cat ", aa_fname[3], " | head"))
aa <- bedtools_merger(bedfile_names =aa_fname , 
                      merge_dist = 1000,
                      col_number=1,
                      o="count",
                      output_filename = "test_output.bed",
                      read_output=T)

#########################################################################################################
#########################################################################################################
# next function starts here


#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         Workflow          ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################


########################################################################################################################
########################################################################################################################
##############################################         ER ChiP         #################################################
########################################################################################################################
########################################################################################################################

# ########################################################################################################################
# # # # # # # # DZIDA_TIMESERIES dataset
# # Collect a set of ER bound regions to be used as the superset of all positive and negative sets
# common_ER_site_list <- list()
# aafiles <- list.files("E_RNA/ESR1/DZIDA_TIMESERIES/")
# for(i in 1:length(aafiles)){
#   common_ER_site_list[[i]] <- readPeakFile(paste0("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/DZIDA_TIMESERIES/", aafiles[i]),as = "data.frame")
# }
# names(common_ER_site_list) <- aafiles
# for(i in 1:length(common_ER_site_list)){
#   print(names(common_ER_site_list)[i])
#   print(nrow(common_ER_site_list[[i]]))
#   print(summary(common_ER_site_list[[i]]$V5 ))
#   print("Number of peaks with > 5 signal/bg ratio")
#   print(sum(common_ER_site_list[[i]]$V5 >= 5))
#   print("#############################")
# }
# 
# # Filter the reads for signal/bg of 5
# common_ER_site_list_filtered_5 <- list()
# for(i in 1:length(common_ER_site_list)){
#   common_ER_site_list_filtered_5[[i]] <- common_ER_site_list[[i]][common_ER_site_list[[i]]$V5 >= 5, ]
# }
# names(common_ER_site_list_filtered_5) <- names(common_ER_site_list)
# # Write bed files for the filtered ones:
# aanames<- unlist(lapply(strsplit(names(common_ER_site_list_filtered_5), split = "\\."), "[[", 1))
# 
# options(scipen=999) # to prevent priniting numbers in scientific format
# for(i in 1:length(common_ER_site_list_filtered_5)){
#   write.table(common_ER_site_list_filtered_5[[i]], file=paste0(aanames[i], "_filtered_5.bed"), quote=F, sep="\t", row.names=F, col.names=F)
# }
# options(scipen=0)
# 
# # excluding some areas? within genes? transcribed regions?
# 
# 
# # check the distribution of the length of peaks:
# common_ER_site_GRanges_list <- list()
# for(i in 1:length(common_ER_site_list)){
#   common_ER_site_GRanges_list[[i]] <- makeGRangesFromDataFrame(df = common_ER_site_list[[i]],
#                                                                keep.extra.columns = T,
#                                                                seqnames.field = "V1", 
#                                                                start.field = "V2",
#                                                                end.field = "V3",
#                                                                ignore.strand = T)
# }
# 
# common_ER_site_GRanges_list_filtered_5 <- list()
# for(i in 1:length(common_ER_site_list)){
#   common_ER_site_GRanges_list_filtered_5[[i]] <- makeGRangesFromDataFrame(df = common_ER_site_list_filtered_5[[i]],
#                                                                keep.extra.columns = T,
#                                                                seqnames.field = "V1", 
#                                                                start.field = "V2",
#                                                                end.field = "V3",
#                                                                ignore.strand = T)
# }
# 
# aa_len_list <- list()
# for(i in 1:length(common_ER_site_list)){
#   aa_len_list[[i]] <- common_ER_site_list[[i]]$V3 - common_ER_site_list[[i]]$V2
# }
# names(aa_len_list) <- names(common_ER_site_list)
# par(mfrow = c(4, 3), mar = c(3,2,4,1))
# for(i in 1:length(aa_len_list)){
#   hist(aa_len_list[[i]],
#        main = paste0(names(aa_len_list)[i]," #" ,length(aa_len_list[[i]])),
#        xlim = c(0, 1800))
# }
# 
# for(i in 1:length(aa_len_list)){
#   print(names(aa_len_list)[i])
#   print(summary(aa_len_list[[i]]))
#   print("#######################")
# }
# # get the union and merge peaks from the first time points: 5, 10, 20, 40, 80
# aa_wd <- getwd()
# setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/DZIDA_TIMESERIES")
# aa_files <- list.files(".")
# 
# system("cat 76108_peaks_ER_E2_5MIN.bed 76107_peaks_ER_E2_10MIN.bed 76106_peaks_ER_E2_20MIN.bed 76105_peaks_ER_E2_40MIN.bed 76104_peaks_ER_E2_80MIN.bed | sortBed -i stdin |mergeBed -i stdin > merged_5_10_20_40_80.bed")
# system("bedtools merge -i merged_5_10_20_40_80.bed -d 100 -c 1 -o count > merged_5_10_20_40_80_100bpmerged.bed")
# system("bedtools merge -i merged_5_10_20_40_80.bed -d 500 -c 1 -o count > merged_5_10_20_40_80_500bpmerged.bed")
# 
# system("cat merged_5_10_20_40_80_100bpmerged.bed | wc -l") # 69879
# 
# system("cat merged_5_10_20_40_80_500bpmerged.bed | wc -l") # 66978
# 
# # system("cat 76108_peaks_ER_E2_5MIN_filtered_5.bed 76107_peaks_ER_E2_10MIN_filtered_5.bed 76106_peaks_ER_E2_20MIN_filtered_5.bed 76105_peaks_ER_E2_40MIN_filtered_5.bed 76104_peaks_ER_E2_80MIN_filtered_5.bed | sortBed -i stdin > sorted_5_10_20_40_80_filtered_5.bed")
# # aa <- readPeakFile("sorted_5_10_20_40_80_filtered_5.bed",as = "data.frame")
# # aa[aa$V1 == "chr11", ]
# 
# 
# system("cat 76108_peaks_ER_E2_5MIN_filtered_5.bed 76107_peaks_ER_E2_10MIN_filtered_5.bed 76106_peaks_ER_E2_20MIN_filtered_5.bed 76105_peaks_ER_E2_40MIN_filtered_5.bed 76104_peaks_ER_E2_80MIN_filtered_5.bed | sortBed -i stdin |mergeBed -i stdin > merged_5_10_20_40_80_filtered_5.bed")
# 
# system("bedtools merge -i merged_5_10_20_40_80_filtered_5.bed -d 100 -c 1 -o count > merged_5_10_20_40_80_100bpmerged_filtered_5.bed")
# system("bedtools merge -i merged_5_10_20_40_80_filtered_5.bed -d 500 -c 1 -o count > merged_5_10_20_40_80_500bpmerged_filtered_5.bed")
# system("bedtools merge -i merged_5_10_20_40_80_filtered_5.bed -d 1000 -c 1 -o count > merged_5_10_20_40_80_1000bpmerged_filtered_5.bed")
# 
# system("cat merged_5_10_20_40_80_100bpmerged_filtered_5.bed | wc -l") # 67652
# 
# system("cat merged_5_10_20_40_80_500bpmerged_filtered_5.bed | wc -l") # 64863
# 
# system("cat merged_5_10_20_40_80_1000bpmerged_filtered_5.bed | wc -l") # 62771
# 
# # look at the length of the peaks:
# aa_100merged <- readPeakFile("merged_5_10_20_40_80_100bpmerged_filtered_5.bed",as = "data.frame")
# aa_500merged <- readPeakFile("merged_5_10_20_40_80_500bpmerged_filtered_5.bed",as = "data.frame")
# aa_1000merged <- readPeakFile("merged_5_10_20_40_80_1000bpmerged_filtered_5.bed",as = "data.frame")
# 
# par(mfrow = c(3, 1), mar= c(4, 4,4,4))
# hist(aa_100merged$V3 - aa_100merged$V2, xlim = c(0,6000), main = paste0("length of 100bp merged. total number: ", nrow(aa_100merged)), freq = F, breaks = 10)
# hist(aa_500merged$V3 - aa_500merged$V2, xlim = c(0, 6000), main = paste0("length of 500bp merged. total number: ", nrow(aa_500merged)), freq = F, breaks = 10)
# hist(aa_1000merged$V3 - aa_1000merged$V2, xlim = c(0, 6000), main = paste0("length of 1000bp merged. total number: ", nrow(aa_1000merged)), freq = F, breaks = 10)
# 
# ###########################
# # Get the union of all ER ChiP peaks as the potential start set (st_er)
# aa_100merged #67652
# 
# ########################################################################################################################
# # # # # # # # Li2013 dataset
# 
# # read ChIP datafrom Li2013 study
# setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/Li2013/")
# aafiles <- list.files()
# ER_ChIP_Li2013 <- list()
# for(i in 1:length(aafiles)){
#   ER_ChIP_Li2013[[i]] <- readPeakFile(aafiles[i],as = "data.frame")
# }
# names(ER_ChIP_Li2013) <- aafiles
# 
# ER_ChIP_Li2013_filtered_5 <- list()
# for(i in 1:length(ER_ChIP_Li2013)){
#   ER_ChIP_Li2013_filtered_5[[i]] <- ER_ChIP_Li2013[[i]][ER_ChIP_Li2013[[i]]$V5 >= 5, ]
# }
# names(ER_ChIP_Li2013_filtered_5) <- names(ER_ChIP_Li2013)
# 
# aanames<- unlist(lapply(strsplit(names(ER_ChIP_Li2013_filtered_5), split = "\\."), "[[", 1))
# options(scipen=999) # to prevent priniting numbers in scientific format
# for(i in 1:length(ER_ChIP_Li2013_filtered_5)){
#   write.table(ER_ChIP_Li2013_filtered_5[[i]], file=paste0(aanames[i],
#                                                           "_filtered_5.bed"), 
#               quote=F, sep="\t", row.names=F, col.names=F)
# }
# options(scipen=0)
# 
# list.files(".")
# system("cat 33100_peaks_e2_1h_rep2_filtered_5.bed 33101_peaks_e2_1h_rep1_filtered_5.bed | sortBed -i stdin |mergeBed -i stdin > merged_ER_e2_1h_union_rep1_rep2_filtered_5.bed")
# system("cat merged_ER_e2_1h_union_rep1_rep2_filtered_5.bed | wc -l") # 34078
# system("bedtools merge -i merged_ER_e2_1h_union_rep1_rep2_filtered_5.bed -d 100 -c 1 -o count > merged_ER_e2_1h_union_rep1_rep2_filtered_5_100bpmerged.bed")
# system("cat merged_ER_e2_1h_union_rep1_rep2_filtered_5_100bpmerged.bed | wc -l") # 34057
# ER_ChIP_Li2013_union_e2_rep1_rep2_100merged <- readPeakFile("merged_ER_e2_1h_union_rep1_rep2_filtered_5_100bpmerged.bed",as = "data.frame")
# 
# system("bedtools intersect -a 33100_peaks_e2_1h_rep2_filtered_5.bed -b 33101_peaks_e2_1h_rep1_filtered_5.bed > intersect_ER_e2_1h_rep1_rep2.bed")
# system("cat 33100_peaks_e2_1h_rep2_filtered_5.bed | wc -l") # 9542
# system("cat 33101_peaks_e2_1h_rep1_filtered_5.bed | wc -l") # 32914
# 
# system("cat intersect_ER_e2_1h_rep1_rep2.bed | wc -l") # 8378
# 
# system("bedtools merge -i intersect_ER_e2_1h_rep1_rep2.bed -d 100 -c 1 -o count > merged_ER_e2_1h_intersect_rep1_rep2_filtered_5_100bpmerged.bed")
# system("cat merged_ER_e2_1h_intersect_rep1_rep2_filtered_5_100bpmerged.bed | wc -l") # 8333
# 
# ## get the intersect with other sets
# aa_a <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/DZIDA_TIMESERIES/merged_5_10_20_40_80_100bpmerged_filtered_5.bed"
# aa_b <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K4ME3/DZIDA_TIMESERIES/merged_H3K4me3_10_20_40_80_filtered_5_100bpmerged.bed"
# aa_c <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/BRUNELLE/merged_H3K27ac_rep1_rep2_Brunelle_filtered_5_100bpmerged.bed"
# aa_d <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K4ME3/DZIDA_TIMESERIES/overlap_ER_100bpmerged_H3K4me3_100bpmerged.bed"
# aa_e <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/overlap_ER_100bpmerged_H3K4me3_100bpmerged_H3K27ac_100bpmerged.bed"
# aa_f <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/overlap_H3K4me3_100bpmerged_H3K27ac_100bpmerged.bed"
# aa_g <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/merged_H2Az_5_10_20_40_80_filtered_5_100bpmerged.bed"
# aa_h <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/Li2013/33099_H3K27AC_CHIPSEQ_E2_LI2013_filtered_5_100bpmerged.bed"
# aa_i <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/union_E21h_rep1_rep2.bed"
# aa_j <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/intersect_E21h_rep1_rep2.bed"
# aa_k <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/non_overlap_E2interes_ETOHinteres.bed"
# aa_l <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/non_overlap_E2inters_ETOHunion.bed"
# aa_m <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/non_overlap_E2union_ETOHunion.bed"
# aa_n <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/non_overlap_E2union_ETOHinters.bed"
# aa_o <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/Li2013/nonoverlap_ER_100bpmerged_Li_H3K27ac_100bpmerged.bed"
# aa_p <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/overlap_ER_100bpmerged_Li_H3k27ac_erna_union_e2.bed"
# aa_q <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/condensin_1248_hg38_active.bed"
# aa_r <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/overlap_ER_100bpmerged_Li_H3k27ac_erna_inters_e2_minus_inters_etoh.bed"
# aa_s <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/Li2013/merged_ER_e2_1h_intersect_rep1_rep2_filtered_5_100bpmerged.bed" # 8333
# aa_t <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/Li2013/merged_ER_e2_1h_union_rep1_rep2_filtered_5_100bpmerged.bed"  # 34057
# aa_u <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/Li2013/overlap_H3k27ac_li_e2_100bpmerged_li_esr1_union.bed"
# 
# 
# system(paste0("bedtools intersect -a ", aa_a, " -b ", aa_t , " > overlap_ER_100bpmerged_li_esr1_union.bed"))
# system("cat overlap_ER_100bpmerged_li_esr1_union.bed | wc -l") #25103
# 
# system(paste0("bedtools intersect -a ", aa_a, " -b ", aa_s , " > overlap_ER_100bpmerged_li_esr1_intersect.bed"))
# system("cat overlap_ER_100bpmerged_li_esr1_intersect.bed | wc -l") #8063
# 
# system(paste0("bedtools intersect -a ", aa_h, " -b ", aa_t , " > overlap_H3k27ac_li_e2_100bpmerged_li_esr1_union.bed"))
# system("cat overlap_H3k27ac_li_e2_100bpmerged_li_esr1_union.bed | wc -l") #3734
# 
# li_H3K27ac_gr <- readPeakFile(aa_h,as = "GRanges")
# li_ER_union_gr <- readPeakFile(aa_t,as = "GRanges")
# li_erna_union_gr <- readPeakFile(aa_i,as = "GRanges")
# li_erna_intersect_gr <- readPeakFile(aa_j,as = "GRanges")
# li_erna_intersect_minus_etoh_intersect_gr <- readPeakFile(aa_k,as = "GRanges")
# 
# vennplot(Sets = list(H3K27ac = li_H3K27ac_gr,
#                      ER_union = li_ER_union_gr
#                      #eRNA_union = li_erna_union_gr
#                      #,eRNA_intersect = li_erna_intersect_gr
#                      ,eRNA_intersect_minus_etoh = li_erna_intersect_minus_etoh_intersect_gr
#                      ),
#          by = "Vennerable")
# 
# Li2013_H3k27ac_esr1_union
# 
# system(paste0("bedtools intersect -a ", aa_u, " -b ", aa_i , " > overlap_H3k27ac_li_e2_100bpmerged_li_esr1_union_erna_union_li.bed"))
# system("cat overlap_H3k27ac_li_e2_100bpmerged_li_esr1_union_erna_union_li.bed | wc -l") #1609
# 
# system(paste0("bedtools intersect -a ", aa_u, " -b ", aa_j , " > overlap_H3k27ac_li_e2_100bpmerged_li_esr1_union_erna_inters_li.bed"))
# system("cat overlap_H3k27ac_li_e2_100bpmerged_li_esr1_union_erna_inters_li.bed | wc -l") #587
# 
# system(paste0("bedtools intersect -a ", aa_u, " -b ", aa_k , " > overlap_H3k27ac_li_e2_100bpmerged_li_esr1_union_erna_inters_minus_etoh_interes.bed"))
# system("cat overlap_H3k27ac_li_e2_100bpmerged_li_esr1_union_erna_inters_minus_etoh_interes.bed | wc -l") #332
# 
# 
# system(paste0("bedtools intersect -a ", aa_h, " -b ", aa_s , " > overlap_H3k27ac_li_e2_100bpmerged_li_esr1_intersect.bed"))
# system("cat overlap_H3k27ac_li_e2_100bpmerged_li_esr1_intersect.bed | wc -l") #1504
# 
# system(paste0("bedtools intersect -a ", aa_t, " -b ", aa_q , " > overlap_li_esr1_union_condensin_active_enhancer.bed"))
# system("cat overlap_li_esr1_union_condensin_active_enhancer.bed | wc -l") #1219
# 
# system(paste0("bedtools intersect -a ", aa_s, " -b ", aa_q , " > overlap_li_esr1_intersect_condensin_active_enhancer.bed"))
# system("cat overlap_li_esr1_intersect_condensin_active_enhancer.bed | wc -l") #1111
# 
# system( paste0("cat ", aa_q, " | wc -l")) # 1247
# 
# system(paste0("bedtools intersect -a ", aa_h, " -b ", aa_q , " > overlap_li_H3k27_condensin_active_enhancer.bed"))
# system("cat overlap_li_H3k27_condensin_active_enhancer.bed | wc -l") #157
# 
# system(paste0("bedtools intersect -a ", aa_c, " -b ", aa_q , " > overlap_brunelle_H3k27_condensin_active_enhancer.bed"))
# system("cat overlap_brunelle_H3k27_condensin_active_enhancer.bed | wc -l") #69
# 
# system(paste0("bedtools intersect -a ", aa_a, " -b ", aa_q , " > overlap_DZIDA_ER_condensin_active_enhancer.bed"))
# system("cat overlap_DZIDA_ER_condensin_active_enhancer.bed | wc -l") #1192
# 
# system(paste0("bedtools intersect -a ", aa_b, " -b ", aa_q , " > overlap_DZIDA_H3K4me3_condensin_active_enhancer.bed"))
# system("cat overlap_DZIDA_H3K4me3_condensin_active_enhancer.bed | wc -l") #171
# 
# system(paste0("bedtools intersect -a ", aa_g, " -b ", aa_q , " > overlap_DZIDA_H2AZ_condensin_active_enhancer.bed"))
# system("cat overlap_DZIDA_H2AZ_condensin_active_enhancer.bed | wc -l") #220
# 
# system(paste0("bedtools intersect -a ", aa_q, " -b ", aa_b , " " ,aa_g, " ", aa_c ," > overlap_condensin_active_enhancer_DZIDA_H3K4me3_H2AZ_brunelle_H3K27ac.bed"))
# system("cat overlap_condensin_active_enhancer_DZIDA_H3K4me3_H2AZ_brunelle_H3K27ac.bed | wc -l") #460
# 
# system(paste0("bedtools intersect -a ", aa_q, " -b ", aa_b , " " ,aa_g, " ", aa_h ," > overlap_condensin_active_enhancer_DZIDA_H3K4me3_H2AZ_li_H3K27ac.bed"))
# system("cat overlap_condensin_active_enhancer_DZIDA_H3K4me3_H2AZ_li_H3K27ac.bed | wc -l") #548
# 
# ########################################################################################################################
# ########################################################################################################################
# ##############################################         H3K4me3         #################################################
# ########################################################################################################################
# ########################################################################################################################
# # Get the union of all H3K4me3 peaks as the potential start set (st_H3k4me3)
# setwd("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K4ME3/DZIDA_TIMESERIES/")
# aafiles <- list.files(".")
# common_H3K4me3_site_list <- list()
# for(i in 1:length(aafiles)){
#   common_H3K4me3_site_list[[i]] <- readPeakFile(aafiles[i],as = "data.frame")
# }
# names(common_H3K4me3_site_list) <- aafiles
# 
# common_H3K4me3_site_list_filtered_5 <- list()
# for(i in 1:length(common_H3K4me3_site_list)){
#   common_H3K4me3_site_list_filtered_5[[i]] <- common_H3K4me3_site_list[[i]][common_H3K4me3_site_list[[i]]$V5 >= 5, ]
# }
# names(common_H3K4me3_site_list_filtered_5) <- names(common_H3K4me3_site_list)
# 
# 
# common_H3K4me3_site_GRanges_list_filtered_5 <- list()
# for(i in 1:length(common_H3K4me3_site_list)){
#   common_H3K4me3_site_GRanges_list_filtered_5[[i]] <- makeGRangesFromDataFrame(df = common_H3K4me3_site_list_filtered_5[[i]],
#                                                                           keep.extra.columns = T,
#                                                                           seqnames.field = "V1", 
#                                                                           start.field = "V2",
#                                                                           end.field = "V3",
#                                                                           ignore.strand = T)
# }
# names(common_H3K4me3_site_GRanges_list_filtered_5) <- names(common_H3K4me3_site_list_filtered_5)
# aanames<- unlist(lapply(strsplit(names(common_H3K4me3_site_GRanges_list_filtered_5), split = "\\."), "[[", 1))
# 
# options(scipen=999) # to prevent priniting numbers in scientific format
# for(i in 1:length(common_H3K4me3_site_list_filtered_5)){
#   write.table(common_H3K4me3_site_list_filtered_5[[i]], file=paste0(aanames[i], "_filtered_5.bed"), quote=F, sep="\t", row.names=F, col.names=F)
# }
# options(scipen=0)
# 
# list.files(".")
# system("cat 76128_peaks_H3K4ME3_E2_10MIN_filtered_5.bed 76127_peaks_H3K4ME3_E2_20MIN_filtered_5.bed 76126_peaks_H3K4ME3_E2_40MIN_filtered_5.bed 76125_peaks_H3K4ME3_E2_80MIN_filtered_5.bed | sortBed -i stdin |mergeBed -i stdin > merged_H3K4me3_10_20_40_80_filtered_5.bed")
# system("cat merged_H3K4me3_10_20_40_80_filtered_5.bed | wc -l") # 40761
# system("bedtools merge -i merged_H3K4me3_10_20_40_80_filtered_5.bed -d 100 -c 1 -o count > merged_H3K4me3_10_20_40_80_filtered_5_100bpmerged.bed")
# system("cat merged_H3K4me3_10_20_40_80_filtered_5_100bpmerged.bed | wc -l") # 37701
# H3K4me3_100merged <- readPeakFile("merged_H3K4me3_10_20_40_80_filtered_5_100bpmerged.bed",as = "data.frame")
# 
# hist(H3K4me3_100merged$V3 - H3K4me3_100merged$V2)
# 
# aa_a <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/DZIDA_TIMESERIES/merged_5_10_20_40_80_100bpmerged_filtered_5.bed"
# aa_b <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K4ME3/DZIDA_TIMESERIES/merged_H3K4me3_10_20_40_80_filtered_5_100bpmerged.bed"
# system(paste0("bedtools intersect -wa -a ", aa_a, " -b ", aa_b , " > overlap_ER_100bpmerged_H3K4me3_100bpmerged.bed"))
# system("cat overlap_ER_100bpmerged_H3K4me3_100bpmerged.bed | wc -l") #9206
# system(paste0("bedtools intersect -v -a ", aa_a, " -b ", aa_b , " > nonoverlap_ER_100bpmerged_H3K4me3_100bpmerged.bed"))
# system("cat nonoverlap_ER_100bpmerged_H3K4me3_100bpmerged.bed | wc -l") #59027
# 
# 
# ########################################################################################################################
# ########################################################################################################################
# ##############################################         H3K27Ac         #################################################
# ########################################################################################################################
# ########################################################################################################################
# #  currently using LI2013 dataset since it's the same study as the li2013 GRO-seq data
# # not using the brunelle data
# # Get the union of all H3k27ac peaks as the potential start set (st_H3K27Ac)
# setwd("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/BRUNELLE")
# aafiles <- list.files(".")
# common_H3k27ac_site_list <- list()
# for(i in 1:length(aafiles)){
#   common_H3k27ac_site_list[[i]] <- readPeakFile(aafiles[i],as = "data.frame")
# }
# names(common_H3k27ac_site_list) <- aafiles
# 
# common_H3k27ac_site_list_filtered_5 <- list()
# for(i in 1:length(common_H3k27ac_site_list)){
#   common_H3k27ac_site_list_filtered_5[[i]] <- common_H3k27ac_site_list[[i]][common_H3k27ac_site_list[[i]]$V5 >= 5, ]
# }
# names(common_H3k27ac_site_list_filtered_5) <- names(common_H3k27ac_site_list)
# 
# aanames<- unlist(lapply(strsplit(names(common_H3k27ac_site_list_filtered_5), split = "\\."), "[[", 1))
# 
# options(scipen=999) # to prevent priniting numbers in scientific format
# for(i in 1:length(common_H3k27ac_site_list_filtered_5)){
#   write.table(common_H3k27ac_site_list_filtered_5[[i]], file=paste0(aanames[i], "_filtered_5.bed"), quote=F, sep="\t", row.names=F, col.names=F)
# }
# options(scipen=0)
# 
# list.files(".")
# system("cat 54486_peaks_H3K27ac_E2_rep1_Brunelle_filtered_5.bed 54486_peaks_H3K27ac_E2_rep2_Brunelle_filtered_5.bed | sortBed -i stdin |mergeBed -i stdin > merged_H3K27ac_rep1_rep2_Brunelle_filtered_5.bed")
# system("cat merged_H3K27ac_rep1_rep2_Brunelle_filtered_5.bed | wc -l") # 22570
# system("bedtools merge -i merged_H3K27ac_rep1_rep2_Brunelle_filtered_5.bed -d 100 -c 1 -o count > merged_H3K27ac_rep1_rep2_Brunelle_filtered_5_100bpmerged.bed")
# system("cat merged_H3K27ac_rep1_rep2_Brunelle_filtered_5_100bpmerged.bed | wc -l") # 21747
# H3K27ac_100merged <- readPeakFile("merged_H3K27ac_rep1_rep2_Brunelle_filtered_5_100bpmerged.bed",as = "data.frame")
# 
# hist(H3K27ac_100merged$V3 - H3K27ac_100merged$V2)
# 
# aa_a <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/DZIDA_TIMESERIES/merged_5_10_20_40_80_100bpmerged_filtered_5.bed"
# aa_b <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K4ME3/DZIDA_TIMESERIES/merged_H3K4me3_10_20_40_80_filtered_5_100bpmerged.bed"
# aa_c <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/BRUNELLE/merged_H3K27ac_rep1_rep2_Brunelle_filtered_5_100bpmerged.bed"
# aa_d <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K4ME3/DZIDA_TIMESERIES/overlap_ER_100bpmerged_H3K4me3_100bpmerged.bed"
# 
# system(paste0("bedtools intersect -wa -a ", aa_a, " -b ", aa_c , " > overlap_ER_100bpmerged_H3K27ac_100bpmerged.bed"))
# system("cat overlap_ER_100bpmerged_H3K27ac_100bpmerged.bed | wc -l") #3518
# system(paste0("bedtools intersect -v -a ", aa_a, " -b ", aa_c , " > nonoverlap_ER_100bpmerged_H3K27ac_100bpmerged.bed"))
# system("cat nonoverlap_ER_100bpmerged_H3K27ac_100bpmerged.bed | wc -l") #64398
# system(paste0("bedtools intersect -a ", aa_d, " -b ", aa_c , " > overlap_ER_100bpmerged_H3K4me3_100bpmerged_H3K27ac_100bpmerged.bed"))
# system("cat overlap_ER_100bpmerged_H3K4me3_100bpmerged_H3K27ac_100bpmerged.bed | wc -l") #2508
# system(paste0("bedtools intersect -a ", aa_b, " -b ", aa_c , " > overlap_H3K4me3_100bpmerged_H3K27ac_100bpmerged.bed"))
# system("cat overlap_H3K4me3_100bpmerged_H3K27ac_100bpmerged.bed | wc -l") #14945
# 
# ########################################################################################################################
# ########################################################################################################################
# ## what about other datasets of H3K27Ac: Li2013 : currently using this dataset
# setwd("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/Li2013/")
# aafiles <- list.files(".")
# 
# aa_li_2013_h3k27 <- list()
# for(i in 1:length(aafiles)){
#   aa_li_2013_h3k27[[i]] <- readPeakFile(aafiles[i],as = "data.frame")
# }
# names(aa_li_2013_h3k27) <- aafiles
# 
# aa_li_2013_h3k27_filtered_5 <- list()
# for(i in 1:length(aa_li_2013_h3k27)){
#   aa_li_2013_h3k27_filtered_5[[i]] <- aa_li_2013_h3k27[[i]][aa_li_2013_h3k27[[i]]$V5 >= 5, ]
# }
# names(aa_li_2013_h3k27_filtered_5) <- names(aa_li_2013_h3k27)
# 
# aanames<- unlist(lapply(strsplit(names(aa_li_2013_h3k27_filtered_5), split = "\\."), "[[", 1))
# 
# options(scipen=999) # to prevent priniting numbers in scientific format
# for(i in 1:length(aa_li_2013_h3k27_filtered_5)){
#   write.table(aa_li_2013_h3k27_filtered_5[[i]], file=paste0(aanames[i], "_filtered_5.bed"), quote=F, sep="\t", row.names=F, col.names=F)
# }
# options(scipen=0)
# list.files(".")
# system("bedtools merge -i 33099_H3K27AC_CHIPSEQ_E2_LI2013_filtered_5.bed -d 100 -c 1 -o count > 33099_H3K27AC_CHIPSEQ_E2_LI2013_filtered_5_100bpmerged.bed")
# aa_li_2013_h3k27_filtered_5_100merged <- readPeakFile("33099_H3K27AC_CHIPSEQ_E2_LI2013_filtered_5_100bpmerged.bed",as = "data.frame")
# H3K27ac_100merged_Li2013 <- aa_li_2013_h3k27_filtered_5_100merged
# hist(aa_li_2013_h3k27_filtered_5_100merged$V3 - aa_li_2013_h3k27_filtered_5_100merged$V2)
# 
# aa_a <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/DZIDA_TIMESERIES/merged_5_10_20_40_80_100bpmerged_filtered_5.bed"
# aa_b <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K4ME3/DZIDA_TIMESERIES/merged_H3K4me3_10_20_40_80_filtered_5_100bpmerged.bed"
# aa_c <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/BRUNELLE/merged_H3K27ac_rep1_rep2_Brunelle_filtered_5_100bpmerged.bed"
# aa_d <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K4ME3/DZIDA_TIMESERIES/overlap_ER_100bpmerged_H3K4me3_100bpmerged.bed"
# aa_e <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/overlap_ER_100bpmerged_H3K4me3_100bpmerged_H3K27ac_100bpmerged.bed"
# aa_f <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/overlap_H3K4me3_100bpmerged_H3K27ac_100bpmerged.bed"
# aa_g <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/merged_H2Az_5_10_20_40_80_filtered_5_100bpmerged.bed"
# aa_h <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/Li2013/33099_H3K27AC_CHIPSEQ_E2_LI2013_filtered_5_100bpmerged.bed"
# 
# system(paste0("bedtools intersect -wa -a ", aa_a, " -b ", aa_h , " > overlap_ER_100bpmerged_Li_H3K27ac_100bpmerged.bed"))
# system("cat overlap_ER_100bpmerged_Li_H3K27ac_100bpmerged.bed | wc -l") #7654
# 
# system(paste0("bedtools intersect -v -a ", aa_a, " -b ", aa_h , " > nonoverlap_ER_100bpmerged_Li_H3K27ac_100bpmerged.bed"))
# system("cat nonoverlap_ER_100bpmerged_Li_H3K27ac_100bpmerged.bed | wc -l") #61045
# 
# system(paste0("bedtools intersect -a ", aa_c, " -b ", aa_h , " > overlap_brunelle_H3k27_Li_H3K27ac_100bpmerged.bed"))
# system("cat overlap_brunelle_H3k27_Li_H3K27ac_100bpmerged.bed | wc -l") #17863
# 
# system(paste0("bedtools intersect -a ", aa_b, " -b ", aa_h , " > overlap_H3K4me3_100bpmerged_Li_H3K27ac_100bpmerged.bed"))
# system("cat overlap_H3K4me3_100bpmerged_Li_H3K27ac_100bpmerged.bed | wc -l") #23164
# 
# system(paste0("bedtools intersect -a ", aa_d, " -b ", aa_h , " > overlap_ER_H3K4me3_100bpmerged_Li_H3K27ac_100bpmerged.bed"))
# system("cat overlap_ER_H3K4me3_100bpmerged_Li_H3K27ac_100bpmerged.bed | wc -l") #5101
# 
# 
# 
# ########################################################################################################################
# ########################################################################################################################
# ##############################################         H2Az         #################################################
# ########################################################################################################################
# ########################################################################################################################
# # Get the union of all H2Az    peaks as the potential start set (st_H2Az)
# setwd("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/")
# aafiles <- list.files(".")
# common_H2Az_site_list <- list()
# for(i in 1:length(aafiles)){
#   common_H2Az_site_list[[i]] <- readPeakFile(aafiles[i],as = "data.frame")
# }
# names(common_H2Az_site_list) <- aafiles
# 
# common_H2Az_site_list_filtered_5 <- list()
# for(i in 1:length(common_H2Az_site_list)){
#   common_H2Az_site_list_filtered_5[[i]] <- common_H2Az_site_list[[i]][common_H2Az_site_list[[i]]$V5 >= 5, ]
# }
# names(common_H2Az_site_list_filtered_5) <- names(common_H2Az_site_list)
# 
# aanames<- unlist(lapply(strsplit(names(common_H2Az_site_list_filtered_5), split = "\\."), "[[", 1))
# 
# options(scipen=999) # to prevent priniting numbers in scientific format
# for(i in 1:length(common_H2Az_site_list_filtered_5)){
#   write.table(common_H2Az_site_list_filtered_5[[i]], file=paste0(aanames[i], "_filtered_5.bed"), quote=F, sep="\t", row.names=F, col.names=F)
# }
# options(scipen=0)
# 
# list.files(".")
# system("cat 76119_peaks_H2AZ_E2_5MIN_filtered_5.bed 76118_peaks_H2AZ_E2_10MIN_filtered_5.bed 76117_peaks_H2AZ_E2_20MIN_filtered_5.bed 76116_peaks_H2AZ_E2_40MIN_filtered_5.bed 76115_peaks_H2AZ_E2_80MIN_filtered_5.bed | sortBed -i stdin |mergeBed -i stdin > merged_H2Az_5_10_20_40_80_filtered_5.bed")
# system("cat merged_H2Az_5_10_20_40_80_filtered_5.bed | wc -l") # 39072
# system("bedtools merge -i merged_H2Az_5_10_20_40_80_filtered_5.bed -d 100 -c 1 -o count > merged_H2Az_5_10_20_40_80_filtered_5_100bpmerged.bed")
# system("cat merged_H2Az_5_10_20_40_80_filtered_5_100bpmerged.bed | wc -l") # 37701
# H2Az_100merged <- readPeakFile("merged_H2Az_5_10_20_40_80_filtered_5_100bpmerged.bed",as = "data.frame")
# 
# hist(H2Az_100merged$V3 - H2Az_100merged$V2)
# 
# aa_a <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/DZIDA_TIMESERIES/merged_5_10_20_40_80_100bpmerged_filtered_5.bed"
# aa_b <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K4ME3/DZIDA_TIMESERIES/merged_H3K4me3_10_20_40_80_filtered_5_100bpmerged.bed"
# aa_c <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/BRUNELLE/merged_H3K27ac_rep1_rep2_Brunelle_filtered_5_100bpmerged.bed"
# aa_d <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K4ME3/DZIDA_TIMESERIES/overlap_ER_100bpmerged_H3K4me3_100bpmerged.bed"
# aa_e <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/overlap_ER_100bpmerged_H3K4me3_100bpmerged_H3K27ac_100bpmerged.bed"
# aa_f <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/overlap_H3K4me3_100bpmerged_H3K27ac_100bpmerged.bed"
# aa_g <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/merged_H2Az_5_10_20_40_80_filtered_5_100bpmerged.bed"
# 
# system(paste0("bedtools intersect -wa -a ", aa_a, " -b ", aa_g , " > overlap_ER_100bpmerged_H2Az_100bpmerged.bed"))
# system("cat overlap_ER_100bpmerged_H2Az_100bpmerged.bed | wc -l") #9380
# system(paste0("bedtools intersect -v -a ", aa_a, " -b ", aa_g , " > nonoverlap_ER_100bpmerged_H2Az_100bpmerged.bed"))
# system("cat nonoverlap_ER_100bpmerged_H2Az_100bpmerged.bed | wc -l") #64398
# system(paste0("bedtools intersect -a ", aa_e, " -b ", aa_g , " > overlap_ER_100bpmerged_H3K4me3_100bpmerged_H3K27ac_100bpmerged_H2Az_100bpmerged.bed"))
# system("cat overlap_ER_100bpmerged_H3K4me3_100bpmerged_H3K27ac_100bpmerged_H2Az_100bpmerged.bed | wc -l") #1884
# system(paste0("bedtools intersect -a ", aa_f, " -b ", aa_g , " > overlap_H3K4me3_100bpmerged_H3K27ac_100bpmerged_H2Az_100bpmerged.bed"))
# system("cat overlap_H3K4me3_100bpmerged_H3K27ac_100bpmerged_H2Az_100bpmerged.bed | wc -l") #14945
# 
# 
# ########################################################################################################################
# ########################################################################################################################
# # look at the length of all mutual intersection, difference and unions
# 
# 
# # Get the intersection of st_er and st_H3k4me3 set_er_p_st_H3k4me3_p
# # Get the difference   of st_er and st_H3k4me3 set_er_p_st_H3k4me3_n
# 
# ########################################################################################################################
# ########################################################################################################################
# ##############################################         GRO-seq         #################################################
# ########################################################################################################################
# ########################################################################################################################
# # get the li2013 erna data
# # read GROSEQ datafrom Li2013 study
# Li2013_GRO_list <- list()
# setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/")
# aafiles <- list.files(".", pattern = "*_hg38.bed")
# for(i in 1:length(aafiles)){
#   Li2013_GRO_list[[i]] <- readPeakFile(aafiles[i], as = "data.frame")
# }
# lapply(Li2013_GRO_list, nrow)
# names(Li2013_GRO_list) <- c("E2_rep1", "E2_rep2", "ETOH_rep1", "ETOH_rep2", "siSMC3_E2_1hr", "siSMC3_EtOH")
# 
# list.files(".")
# system("cat SRR816998-1_divergent_classifications.bednoinvalidlines_E2_1h_hg38.bed SRR816999-1_divergent_classifications.bednoinvalidlines_E21h_hg38.bed | sortBed -i stdin |mergeBed -i stdin > union_E21h_rep1_rep2.bed")
# system("cat union_E21h_rep1_rep2.bed | wc -l") # 35153
# aa_union_erna_e2_rep1_rep2 <- readPeakFile("union_E21h_rep1_rep2.bed",as = "data.frame")
# 
# system("cat SRR817000-1_divergent_classifications.bednoinvalidlines_ETOH_1h_hg38.bed SRR817001-1_divergent_classifications.bednoinvalidlines_ETOH_1h_hg38.bed | sortBed -i stdin |mergeBed -i stdin > union_ETOH1h_rep1_rep2.bed")
# system("cat union_ETOH1h_rep1_rep2.bed | wc -l") # 38769
# aa_union_erna_etoh_rep1_rep2 <- readPeakFile("union_ETOH1h_rep1_rep2.bed",as = "data.frame")
# 
# #### do I need to exclude known transcripts as mentioned in:  A Pan-Cancer Analysis of Enhancer Expression in Nearly 9000 Patient Samples
# sum(aa_union_erna_e2_rep1_rep2$V3 - aa_union_erna_e2_rep1_rep2$V2 > 2000)
# summary(Li2013_GRO_list$E2_rep1$V3 - Li2013_GRO_list$E2_rep1$V2 > 2000)
# summary(Li2013_GRO_list$E2_rep2$V3 - Li2013_GRO_list$E2_rep2$V2 > 2000)
# 
# system("bedtools intersect -a SRR816998-1_divergent_classifications.bednoinvalidlines_E2_1h_hg38.bed -b SRR816999-1_divergent_classifications.bednoinvalidlines_E21h_hg38.bed > intersect_E21h_rep1_rep2.bed")
# system("cat intersect_E21h_rep1_rep2.bed | wc -l") # 11524
# aa_intersect_erna_e2_rep1_rep2 <- readPeakFile("intersect_E21h_rep1_rep2.bed",as = "data.frame")
# summary(aa_intersect_erna_e2_rep1_rep2$V3 - aa_intersect_erna_e2_rep1_rep2$V2)
# 
# system("bedtools intersect -a SRR817000-1_divergent_classifications.bednoinvalidlines_ETOH_1h_hg38.bed -b SRR817001-1_divergent_classifications.bednoinvalidlines_ETOH_1h_hg38.bed > intersect_ETOH1h_rep1_rep2.bed")
# system("cat intersect_ETOH1h_rep1_rep2.bed | wc -l") # 11646
# aa_intersect_erna_etoh_rep1_rep2 <- readPeakFile("intersect_ETOH1h_rep1_rep2.bed",as = "data.frame")
# summary(aa_intersect_erna_etoh_rep1_rep2$V3 - aa_intersect_erna_etoh_rep1_rep2$V2)
# 
# 
# system("bedtools intersect -v -wa -a intersect_E21h_rep1_rep2.bed -b intersect_ETOH1h_rep1_rep2.bed > non_overlap_E2interes_ETOHinteres.bed")
# system("cat non_overlap_E2interes_ETOHinteres.bed | wc -l") # 4941
# 
# system("bedtools intersect -v  -a union_E21h_rep1_rep2.bed -b union_ETOH1h_rep1_rep2.bed > non_overlap_E2union_ETOHunion.bed")
# system("cat non_overlap_E2union_ETOHunion.bed | wc -l") # 13129
# 
# system("bedtools intersect -v  -a intersect_E21h_rep1_rep2.bed -b union_ETOH1h_rep1_rep2.bed > non_overlap_E2inters_ETOHunion.bed")
# system("cat non_overlap_E2inters_ETOHunion.bed | wc -l") # 1531
# 
# system("bedtools intersect -v  -a union_E21h_rep1_rep2.bed -b intersect_ETOH1h_rep1_rep2.bed > non_overlap_E2union_ETOHinters.bed")
# system("cat non_overlap_E2union_ETOHinters.bed | wc -l") # 25757
# 
# 
# aa_a <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/DZIDA_TIMESERIES/merged_5_10_20_40_80_100bpmerged_filtered_5.bed"
# aa_b <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K4ME3/DZIDA_TIMESERIES/merged_H3K4me3_10_20_40_80_filtered_5_100bpmerged.bed"
# aa_c <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/BRUNELLE/merged_H3K27ac_rep1_rep2_Brunelle_filtered_5_100bpmerged.bed"
# aa_d <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K4ME3/DZIDA_TIMESERIES/overlap_ER_100bpmerged_H3K4me3_100bpmerged.bed"
# aa_e <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/overlap_ER_100bpmerged_H3K4me3_100bpmerged_H3K27ac_100bpmerged.bed"
# aa_f <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/overlap_H3K4me3_100bpmerged_H3K27ac_100bpmerged.bed"
# aa_g <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H2AZ/DZIDA_TIMESERIES/merged_H2Az_5_10_20_40_80_filtered_5_100bpmerged.bed"
# aa_h <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/Li2013/33099_H3K27AC_CHIPSEQ_E2_LI2013_filtered_5_100bpmerged.bed"
# aa_i <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/union_E21h_rep1_rep2.bed"
# aa_j <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/intersect_E21h_rep1_rep2.bed"
# aa_k <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/non_overlap_E2interes_ETOHinteres.bed"
# aa_l <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/non_overlap_E2inters_ETOHunion.bed"
# aa_m <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/non_overlap_E2union_ETOHunion.bed"
# aa_n <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/non_overlap_E2union_ETOHinters.bed"
# aa_o <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/Li2013/nonoverlap_ER_100bpmerged_Li_H3K27ac_100bpmerged.bed"
# aa_p <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/overlap_ER_100bpmerged_Li_H3k27ac_erna_union_e2.bed"
# aa_q <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/condensin_1248_hg38_active.bed"
# aa_r <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/overlap_ER_100bpmerged_Li_H3k27ac_erna_inters_e2_minus_inters_etoh.bed"
# system(paste0("bedtools intersect -wa -a ", aa_a, " -b ", aa_i , " > overlap_ER_100bpmerged_erna_union_e2.bed"))
# system("cat overlap_ER_100bpmerged_erna_union_e2.bed | wc -l") #6249
# 
# system(paste0("bedtools intersect -wa -a ", aa_a, " -b ", aa_j , " > overlap_ER_100bpmerged_erna_intersection_e2.bed"))
# system("cat overlap_ER_100bpmerged_erna_intersection_e2.bed | wc -l") #2197
# 
# system(paste0("bedtools intersect -wa -a ", aa_a, " -b ", aa_m , " > overlap_ER_100bpmerged_erna_union_e2_minus_union_etoh.bed"))
# system("cat overlap_ER_100bpmerged_erna_union_e2_minus_union_etoh.bed | wc -l") #2192
# 
# system(paste0("bedtools intersect -wa -a ", aa_a, " -b ", aa_n , " > overlap_ER_100bpmerged_erna_union_e2_minus_inters_etoh.bed"))
# system("cat overlap_ER_100bpmerged_erna_union_e2_minus_inters_etoh.bed | wc -l") #4239
# 
# system(paste0("bedtools intersect -wa -a ", aa_a, " -b ", aa_k , " > overlap_ER_100bpmerged_erna_inters_e2_minus_inters_etoh.bed"))
# system("cat overlap_ER_100bpmerged_erna_inters_e2_minus_inters_etoh.bed | wc -l") #1203
# 
# system(paste0("bedtools intersect -wa -a ", aa_o, " -b ", aa_i , " > overlap_ER_100bpmerged_Li_H3k27ac_erna_union_e2.bed"))
# system("cat overlap_ER_100bpmerged_Li_H3k27ac_erna_union_e2.bed | wc -l") #3300
# 
# system(paste0("bedtools intersect -wa -a ", aa_o, " -b ", aa_k , " > overlap_ER_100bpmerged_Li_H3k27ac_erna_inters_e2_minus_inters_etoh.bed"))
# system("cat overlap_ER_100bpmerged_Li_H3k27ac_erna_inters_e2_minus_inters_etoh.bed | wc -l") #625
# 
# 
# # compare the 1248 enhancers from condensin 1 2 complex paper with the obtained 1203 
# aa_condensin_1248 <- read.delim("Condensin_paper_1248_Active_enh.csv", sep = ",", header = F)
# rownames(aa_condensin_1248) <- aa_condensin_1248[, 1]
# aa_condensin_1248[, 1] <- NULL
# write.table(aa_condensin_1248, file="aa_condensin_1248.bed", quote=F, sep="\t", row.names=F, col.names=F)
# # converted to hg38:
# system(paste0("bedtools intersect -wa -a ", aa_p, " -b ", aa_q , " > overlap_ER_100bpmerged_Li_H3k27ac_erna_union_e2_condensin.bed"))
# system("cat overlap_ER_100bpmerged_Li_H3k27ac_erna_union_e2_condensin.bed | wc -l") #423
# 
# system(paste0("bedtools intersect -wa -a ", aa_r, " -b ", aa_q , " > overlap_ER_100bpmerged_Li_H3k27ac_erna_inters_e2_minus_interes_etoh_condensin.bed"))
# system("cat overlap_ER_100bpmerged_Li_H3k27ac_erna_inters_e2_minus_interes_etoh_condensin.bed | wc -l") #222
# 
# system(paste0("bedtools intersect -a ", aa_i, " -b ", aa_q , " > overlap_erna_union_e2_condensin.bed"))
# system("cat overlap_erna_union_e2_condensin.bed | wc -l") #789
# 
# system(paste0("bedtools intersect -a ", aa_j, " -b ", aa_q , " > overlap_erna_inters_e2_condensin.bed"))
# system("cat overlap_erna_inters_e2_condensin.bed | wc -l") #492
# 

########################################################################################################################

# Decided to work with datasets based on H3K27Ac, e-RNA and ER ChIP from Li 2013 study for now
# FUNCTION THAT CAN BE USED MAYBE WITH MODIFICATION
WriteFastaOfBag
AddpsudoRow
RemoveDegenRow
PWMtoCount
MotifWriter
ExpressionWriter
LogFoldTransform
Hal_job_writer
EnhancerChopper_rangeMat
EnhancerChopper
EnhancerChopper_wrapper
RangeChecker
EnsembleParConstructor_Multi_enh
FactorInfoWriter
par_ff_lb_ub_coop_Writer_multi_enh
initialPar_Writer_multi_enh
freefixWriter_multi_enh
upper_lower_bound_Writer_multi_enh
coop_Writer_Multi_enh
Hal_job_writer_Multi_enh
################################################################################################################
li_H3K27ac <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/Li2013/33099_H3K27AC_CHIPSEQ_E2_LI2013_filtered_5_100bpmerged.bed"
li_ER_union <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/Li2013/merged_ER_e2_1h_union_rep1_rep2_filtered_5_100bpmerged.bed"  # 34057
li_erna_intersect <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/e_RNA_rep1_rep2_intersecting_reads.bed"

# getting the e-RNAs that appeared in each of the reps that overlap with the other rep
aafnames <- list.files("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/" , pattern = "*E2_1h_hg38.bed", full.names = T)
aa_tmp_erna_1_2_inters <- bedtools_intersect(bedfile_names = aafnames, wa=F, wb=F, loj=F, wo=F, wao=F,
                                             u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                             S=F, sorted=F, merge_after_distance=100, 
                                             output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/tmp_E2_1_2_intersect.bed",
                                             read_output = T)

aa_tmp_erna_2_1_inters <- bedtools_intersect(bedfile_names = rev(aafnames), wa=F, wb=F, loj=F, wo=F, wao=F,
                                             u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                             S=F, sorted=F, merge_after_distance=100, 
                                             output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/tmp_E2_2_1_intersect.bed",
                                             read_output = T)

aa_fname <- c("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/tmp_E2_1_2_intersect.bed", 
              "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/tmp_E2_2_1_intersect.bed")

aa <- bedtools_merger(bedfile_names =aa_fname , 
                                     merge_dist = 100,
                                     col_number=1,
                                     o="count",
                                     output_filename = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/e_RNA_rep1_rep2_intersecting_reads.bed",
                                     read_output=T)
li_erna_intersect <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/e_RNA_rep1_rep2_intersecting_reads.bed"
aa_tmp_erna_1_2_union <- bedtools_merger(bedfile_names = c("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/SRR816998-1_divergent_classifications.bednoinvalidlines_E2_1h_hg38.bed",
                                                              "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/SRR816999-1_divergent_classifications.bednoinvalidlines_E2_1h_hg38.bed"),
                                         merge_dist = 100,  
                                            output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/2/tmp_erna_1_2_union.bed",
                                         read_output = T)

li_erna_union <-"~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/2/tmp_erna_1_2_union.bed"

#### take intersect of the ( intersect of the two eRNA replicates : li_erna_intersect), (union of the two ER chip: li_ER_union), (H3k27ac : li_ER_union) as the positive set
Positive_set_gr_list <- list()
aa_tmp_H3K27ac_ER <- bedtools_intersect(bedfile_names = c(li_ER_union, li_H3K27ac), wa=F, wb=F, loj=F, wo=F, wao=F,
                                        u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                        S=F, sorted=F, merge_after_distance=100, 
                                        output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/1/tmp_H3K27ac_ER.bed", read_output = T)
aa_tmp_H3K27ac_ER_name <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/1/tmp_H3K27ac_ER.bed"
Positive_set_gr_list[[1]] <- bedtools_intersect(bedfile_names = c(aa_tmp_H3K27ac_ER_name,  li_erna_intersect), wa=F, wb=F, loj=F, wo=F, wao=F,
                                                u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                                S=F, sorted=F, merge_after_distance=100, 
                                                output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/1/positive_1.bed", read_output = T)
# take difference of (intersect of (union of the two ER chip) and H3k27ac ) and ( intersect of the two eRNA replicates) as negative set
Negative_set_gr_list <- list()
Negative_set_gr_list[[1]] <- bedtools_intersect(bedfile_names = c(aa_tmp_H3K27ac_ER_name,  li_erna_intersect), wa=F, wb=F, loj=F, wo=F, wao=F,
                                                u=F, c=F, v=T, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                                S=F, sorted=F, merge_after_distance=100, 
                                                output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/1/negative_1.bed", read_output = T)

Accepted_chromosomes <- paste0("chr", c(c(1:22), "X", "Y"))
aadiff <- setdiff(levels(seqnames(Negative_set_gr_list[[1]])) , Accepted_chromosomes)
Negative_set_gr_list[[1]] <- Negative_set_gr_list[[1]][! seqnames(Negative_set_gr_list[[1]]) %in% aadiff]
options(scipen=999)
write.table(Negative_set_gr_list[[1]], 
              file="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/1/negative_1.bed",
              quote=F, sep="\t", row.names=F, col.names=F)

options(scipen=0)
li_H3K27ac_gr <- readPeakFile(li_H3K27ac, as = "GRanges")
li_ER_union_gr <- readPeakFile(li_ER_union, as = "GRanges")
li_erna_intersect_gr <- readPeakFile(li_erna_intersect, as = "GRanges")
vennplot(Sets = list(H3K27ac = li_H3K27ac_gr,
                     ER_union = li_ER_union_gr
                     #eRNA_union = li_erna_union_gr
                     ,eRNA_intersect = li_erna_intersect_gr
                     #,eRNA_intersect_minus_etoh = li_erna_intersect_minus_etoh_intersect_gr
                     ),
         by = "Vennerable")
################################################################################################################
# create a positive_negative_Set number 2: positive is the same as one, negative is intersect( union(ER), H3k27) - union(eRNA reps)
Positive_set_gr_list[[2]] <- bedtools_intersect(bedfile_names = c(aa_tmp_H3K27ac_ER_name,  li_erna_intersect), wa=F, wb=F, loj=F, wo=F, wao=F,
                                                u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                                S=F, sorted=F, merge_after_distance=100, 
                                                output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/2/positive_2.bed", read_output = T)
Negative_set_gr_list[[2]] <- bedtools_intersect(bedfile_names = c(aa_tmp_H3K27ac_ER_name,  li_erna_union), wa=F, wb=F, loj=F, wo=F, wao=F,
                                                u=F, c=F, v=T, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                                S=F, sorted=F, merge_after_distance=100, 
                                                output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/2/negative_2.bed", read_output = T)
Accepted_chromosomes <- paste0("chr", c(c(1:22), "X", "Y"))
aadiff <- setdiff(levels(seqnames(Negative_set_gr_list[[2]])) , Accepted_chromosomes)
Negative_set_gr_list[[2]] <- Negative_set_gr_list[[2]][! seqnames(Negative_set_gr_list[[2]]) %in% aadiff]
options(scipen=999)
write.table(Negative_set_gr_list[[2]], 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/2/negative_2.bed",
            quote=F, sep="\t", row.names=F, col.names=F)

options(scipen=0)
################################################################################################################
# create a positive_negative_Set number 3: 
aa_tmp_H3K27ac_ER <- bedtools_intersect(bedfile_names = c(li_ER_union, li_H3K27ac), wa=F, wb=F, loj=F, wo=F, wao=F,
                                        u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                        S=F, sorted=F, merge_after_distance=100, 
                                        output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/3/tmp_H3K27ac_ER.bed", read_output = T)
aa_tmp_H3K27ac_ER_name <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/3/tmp_H3K27ac_ER.bed"
aa_tmp_ER_no_H3k27 <- bedtools_intersect(bedfile_names = c(li_ER_union, li_H3K27ac), wa=F, wb=F, loj=F, wo=F, wao=F,
                                        u=F, c=F, v=T, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                        S=F, sorted=F, merge_after_distance=100, 
                                        output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/3/tmp_ER_no_H3k27.bed", read_output = T)
aa_tmp_ER_no_H3k27_name <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/3/tmp_ER_no_H3k27.bed"
Positive_set_gr_list[[3]] <- bedtools_intersect(bedfile_names = c(aa_tmp_H3K27ac_ER_name,  li_erna_intersect), wa=F, wb=F, loj=F, wo=F, wao=F,
                                                u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                                S=F, sorted=F, merge_after_distance=100, 
                                                output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/3/positive_3.bed", read_output = T)
Negative_set_gr_list[[3]] <- bedtools_intersect(bedfile_names = c(aa_tmp_ER_no_H3k27_name,  li_erna_union), wa=F, wb=F, loj=F, wo=F, wao=F,
                                                u=F, c=F, v=T, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                                S=F, sorted=F, merge_after_distance=100, 
                                                output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/3/negative_3.bed", read_output = T)
Accepted_chromosomes <- paste0("chr", c(c(1:22), "X", "Y"))
aadiff <- setdiff(levels(seqnames(Negative_set_gr_list[[3]])) , Accepted_chromosomes)
Negative_set_gr_list[[3]] <- Negative_set_gr_list[[3]][! seqnames(Negative_set_gr_list[[3]]) %in% aadiff]
options(scipen=999)
write.table(Negative_set_gr_list[[3]], 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/3/negative_3.bed",
            quote=F, sep="\t", row.names=F, col.names=F)

options(scipen=0)
################################################################################################################
# create a positive_negative_Set number 4: Positive set: using the e-RNA peaks that only appeared after e2 as the positive set, removing the ones that were also positive before e2 stimulation (intersect of two ETOH reps) 
# Negative set is the same as negative set 3: the ones with ER and H3k27 signal that don't have eRNA produced in either replicates of E2 treatment.
aa_tmp_H3K27ac_ER_name <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/1/tmp_H3K27ac_ER.bed"
li_erna_intersect <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/e_RNA_rep1_rep2_intersecting_reads.bed"

"union_ETOH1h_rep1_rep2.bed"
aa_non_E2_eRNA <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/intersect_ETOH1h_rep1_rep2.bed"
aa_pure_e2_eRNA <- bedtools_intersect(bedfile_names = c(li_erna_intersect,  aa_non_E2_eRNA),
                                      wa=F, wb=F, loj=F, wo=F, wao=F,
                                      u=F, c=F, v=T, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                      S=F, sorted=F, merge_after_distance=100, 
                                      output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/4/tmp_instersect_eRNA_OnlyAfterE2.bed",
                                      read_output = T)

aa_tmp_erna_only_afterE2 <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/4/tmp_instersect_eRNA_OnlyAfterE2.bed"
aatest <- bedtools_intersect(bedfile_names = c(aa_tmp_H3K27ac_ER_name,  aa_tmp_erna_only_afterE2),
                             wa=F, wb=F, loj=F, wo=F, wao=F,
                             u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                             S=F, sorted=F, merge_after_distance=100, 
                             output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/4/positive_4.bed", read_output = T)
aatest2 <- bedtools_intersect(bedfile_names = c(aa_tmp_H3K27ac_ER_name,  li_erna_intersect), wa=F, wb=F, loj=F, wo=F, wao=F,
                              u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                              S=F, sorted=F, merge_after_distance=100, 
                              output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/4/positive_4.bed", read_output = T)

Positive_set_gr_list[[4]] <- bedtools_intersect(bedfile_names = c(aa_tmp_H3K27ac_ER_name,  aa_tmp_erna_only_afterE2), 
                                                wa=F, wb=F, loj=F, wo=F, wao=F,
                                                u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                                S=F, sorted=F, merge_after_distance=100, 
                                                output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/4/positive_4.bed", read_output = T)
Negative_set_gr_list[[4]] <- bedtools_intersect(bedfile_names = c(aa_tmp_H3K27ac_ER_name,  li_erna_union),
                                                wa=F, wb=F, loj=F, wo=F, wao=F,
                                                u=F, c=F, v=T, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                                S=F, sorted=F, merge_after_distance=100, 
                                                output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/4/negative_4.bed", read_output = T)
Accepted_chromosomes <- paste0("chr", c(c(1:22), "X", "Y"))
aadiff <- setdiff(levels(seqnames(Negative_set_gr_list[[4]])) , Accepted_chromosomes)
Negative_set_gr_list[[4]] <- Negative_set_gr_list[[4]][! seqnames(Negative_set_gr_list[[4]]) %in% aadiff]
options(scipen=999)
write.table(Negative_set_gr_list[[4]], 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/4/negative_4.bed",
            quote=F, sep="\t", row.names=F, col.names=F)

options(scipen=0)


li_H3K27ac_gr <- readPeakFile(li_H3K27ac, as = "GRanges")
li_ER_union_gr <- readPeakFile(li_ER_union, as = "GRanges")
li_erna_only_gr <- readPeakFile(aa_tmp_erna_only_afterE2, as = "GRanges")
vennplot(Sets = list(H3K27Ac = li_H3K27ac_gr,
                     ER = li_ER_union_gr
                     #eRNA_union = li_erna_union_gr
                     ,eRNA = li_erna_only_gr
                     #,eRNA_intersect_minus_etoh = li_erna_intersect_minus_etoh_intersect_gr
),
by = "Vennerable")

## testing to find the intensity field
aatst <- read.table("~/Downloads/SRR816998-1_divergent_classifications.bednoinvalidlines.bed",
                    stringsAsFactors = F, header = F)
aaa <- aatst$V4
aaa <- strsplit(aaa, split = "\\|")
aaa <- lapply(aaa, "[[", 2)
aaa2 <- lapply(aaa, strsplit, split= ",")
aaa3 <- lapply(aaa2, unlist)
aaa4 <- do.call(rbind, aaa3)
aaa5 <- apply(aaa4, 1, as.numeric)
aaa5 <- t(aaa5)
summary(aaa5)
aaziad <- which(aaa5[, 2] > 500)
aakam <- which(aaa5[, 2] < 25)
aatst[aaziad[1],]
aatst[aakam[2],]
aaa5[aaziad[1], ]
################################################################################################################
# obtain the sequence for positive and negative sets
genome38 <- BSgenome.Hsapiens.UCSC.hg38
knownGene_txdb_38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
Positive_set_seq_list <- list()
Negative_set_seq_list <- list()

Positive_set_seq_list[[1]] <- getSeq(x = genome38, names = Positive_set_gr_list[[1]])
Negative_set_seq_list[[1]] <- getSeq(x = genome38, names = Negative_set_gr_list[[1]])

par(mfrow = c(2, 1), mar = c(4,4,4,4))
hist(width(Positive_set_seq_list[[1]]),main = "positive set length")
hist(width(Negative_set_seq_list[[1]]), main = "negative set length")

# create a list where all input sequences are 1000 bp long -- if shorter add flank, if longer cut flank
Positive_set_seq_list_1000bp <- list()
Negative_set_seq_list_1000bp <- list()
Positive_set_gr_list_1000bp <- list()
Negative_set_gr_list_1000bp <- list()

aa_Cur_width <- width(Positive_set_gr_list[[1]])
aa_change <- (1000 - aa_Cur_width)/2
aaranges <- Positive_set_gr_list[[1]]
ranges(aaranges) <- ranges(aaranges) + aa_change
Positive_set_gr_list_1000bp[[1]] <- aaranges

aa_Cur_width <- width(Negative_set_gr_list[[1]])
aa_change <- (1000 - aa_Cur_width)/2
aaranges <- Negative_set_gr_list[[1]]
ranges(aaranges) <- ranges(aaranges) + aa_change
Negative_set_gr_list_1000bp[[1]] <- aaranges

Positive_set_seq_list_1000bp[[1]] <- getSeq(x = genome38, names = Positive_set_gr_list_1000bp[[1]])
Negative_set_seq_list_1000bp[[1]] <- getSeq(x = genome38, names = Negative_set_gr_list_1000bp[[1]])
save(list = c("Positive_set_seq_list_1000bp", "Negative_set_seq_list_1000bp"), file = "Positive_Negative_1000bp_1.RData")

# creating 1000bp versions of the second positive_negative_set
aa_Cur_width <- width(Positive_set_gr_list[[2]])
aa_change <- (1000 - aa_Cur_width)/2
aaranges <- Positive_set_gr_list[[2]]
ranges(aaranges) <- ranges(aaranges) + aa_change
Positive_set_gr_list_1000bp[[2]] <- aaranges

aa_Cur_width <- width(Negative_set_gr_list[[2]])
aa_change <- (1000 - aa_Cur_width)/2
aaranges <- Negative_set_gr_list[[2]]
ranges(aaranges) <- ranges(aaranges) + aa_change
Negative_set_gr_list_1000bp[[2]] <- aaranges

Positive_set_seq_list_1000bp[[2]] <- getSeq(x = genome38, names = Positive_set_gr_list_1000bp[[2]])
Negative_set_seq_list_1000bp[[2]] <- getSeq(x = genome38, names = Negative_set_gr_list_1000bp[[2]])

aamatch <- match(Negative_set_seq_list_1000bp[[2]], Negative_set_seq_list_1000bp[[1]])
################################################################################################################
# sequence for thirs set
aa_Cur_width <- width(Positive_set_gr_list[[3]])
aa_change <- (1000 - aa_Cur_width)/2
aaranges <- Positive_set_gr_list[[3]]
ranges(aaranges) <- ranges(aaranges) + aa_change
Positive_set_gr_list_1000bp[[3]] <- aaranges

aa_Cur_width <- width(Negative_set_gr_list[[3]])
aa_change <- (1000 - aa_Cur_width)/2
aaranges <- Negative_set_gr_list[[3]]
ranges(aaranges) <- ranges(aaranges) + aa_change
Negative_set_gr_list_1000bp[[3]] <- aaranges

Positive_set_seq_list_1000bp[[3]] <- getSeq(x = genome38, names = Positive_set_gr_list_1000bp[[3]])
Negative_set_seq_list_1000bp[[3]] <- getSeq(x = genome38, names = Negative_set_gr_list_1000bp[[3]])
################################################################################################################
# sequence for 4th set
aa_Cur_width <- width(Positive_set_gr_list[[4]])
aa_change <- (1000 - aa_Cur_width)/2
aaranges <- Positive_set_gr_list[[4]]
ranges(aaranges) <- ranges(aaranges) + aa_change
Positive_set_gr_list_1000bp[[4]] <- aaranges

aa_Cur_width <- width(Negative_set_gr_list[[4]])
aa_change <- (1000 - aa_Cur_width)/2
aaranges <- Negative_set_gr_list[[4]]
ranges(aaranges) <- ranges(aaranges) + aa_change
Negative_set_gr_list_1000bp[[4]] <- aaranges

Positive_set_seq_list_1000bp[[4]] <- getSeq(x = genome38, names = Positive_set_gr_list_1000bp[[4]])
Negative_set_seq_list_1000bp[[4]] <- getSeq(x = genome38, names = Negative_set_gr_list_1000bp[[4]])

length(intersect(Positive_set_seq_list_1000bp[[4]], Positive_set_seq_list_1000bp[[2]]))
length(intersect(Negative_set_seq_list_1000bp[[4]], Negative_set_seq_list_1000bp[[2]]))

################################################################################################################
# scanning the sequences for motifs
Positive_set_motifscore <- list()
Positive_set_motifscore[[1]] <- list()
Negative_set_motifscore <- list()
Negative_set_motifscore[[1]] <- list()
#motif list
#TF.motifs.Expanded
TF.motifs.Expanded_t <- lapply(TF.motifs.Expanded, t)
TF.motifs.Expanded_pseudo <- lapply(TF.motifs.Expanded, AddpsudoRow)
TF.motifs.Expanded_pseudo_t <- lapply(TF.motifs.Expanded_pseudo, t)

# for(i in 1:length(Positive_set_seq_list[[1]])){
#   print(paste("positive seq nu ", i, " out of ", length(Positive_set_seq_list[[1]])))
#   Positive_set_motifscore[[1]][[i]] <- MotifScore(seq = as.character(Positive_set_seq_list[[1]][i]),
#                                                   bg = c(0.25, 0.25, 0.25, 0.25),
#                                                   motifList = TF.motifs.Expanded_pseudo_t)
# }
# for(i in 1:length(Negative_set_seq_list[[1]])){
#   print(paste("negative seq nu ", i, " out of ", length(Negative_set_seq_list[[1]])))
#   Negative_set_motifscore[[1]][[i]] <- MotifScore(seq = as.character(Negative_set_seq_list[[1]][i]),
#                                                   bg = c(0.25, 0.25, 0.25, 0.25),
#                                                   motifList = TF.motifs.Expanded_pseudo_t)
# 
# }
# load("motifScore_psuedo.RData")
names(Positive_set_motifscore[[1]]) <- paste0("pos_", c(1:length(Positive_set_motifscore[[1]])))
names(Negative_set_motifscore[[1]]) <- paste0("neg_", c(1:length(Negative_set_motifscore[[1]])))
################################################################################################################
# computing motif scores for the 1000 bp sequences
Positive_set_motifscore_1000bp <- list()
Positive_set_motifscore_1000bp[[1]] <- list()
Negative_set_motifscore_1000bp <- list()
Negative_set_motifscore_1000bp[[1]] <- list()

# for(i in 1:length(Positive_set_seq_list_1000bp[[1]])){
#   print(paste("positive seq nu ", i, " out of ", length(Positive_set_seq_list_1000bp[[1]])))
#   Positive_set_motifscore_1000bp[[1]][[i]] <- MotifScore(seq = as.character(Positive_set_seq_list_1000bp[[1]][i]),
#                                                   bg = c(0.25, 0.25, 0.25, 0.25),
#                                                   motifList = TF.motifs.Expanded_pseudo_t)
# }
# for(i in 1:length(Negative_set_seq_list_1000bp[[1]])){
#   print(paste("negative seq nu ", i, " out of ", length(Negative_set_seq_list_1000bp[[1]])))
#   Negative_set_motifscore_1000bp[[1]][[i]] <- MotifScore(seq = as.character(Negative_set_seq_list_1000bp[[1]][i]),
#                                                   bg = c(0.25, 0.25, 0.25, 0.25),
#                                                   motifList = TF.motifs.Expanded_pseudo_t)
# 
# }
# read half from hal, half from saved RData
#load("MotifScore_Results_Positive_full_Negative_1_574_1000bp.RData")
Positive_set_seq_list_char_1000bp <- list()
Positive_set_seq_list_char_1000bp[[1]] <- as.character(Positive_set_seq_list_1000bp[[1]])
names(Positive_set_seq_list_char_1000bp[[1]]) <- paste("pos", c(1:length(Positive_set_seq_list_1000bp[[1]])), sep = "_")

# read the rest of negatives from hal output
# for(i in 1027:length(Negative_set_seq_list_1000bp[[1]])){
#   print(i)
#   load(paste0("Negative_Set_hal/MotifScore_neg_", i, ".RData"))
#   Negative_set_motifscore_1000bp[[1]][[i]] <- my_motif_score
#   rm("my_motif_score")
# }

Negative_set_seq_list_char_1000bp <- list()
Negative_set_seq_list_char_1000bp[[1]] <- as.character(Negative_set_seq_list_1000bp[[1]])
names(Negative_set_seq_list_char_1000bp[[1]]) <- paste("neg", c(1:length(Negative_set_seq_list_1000bp[[1]])), sep = "_")

Negative_set_seq_list_char_1000bp[[2]] <- as.character(Negative_set_seq_list_1000bp[[2]])
aamatch <- match(Negative_set_seq_list_1000bp[[2]], Negative_set_seq_list_1000bp[[1]])
names(Negative_set_seq_list_char_1000bp[[2]]) <- names(Negative_set_seq_list_char_1000bp[[1]])[aamatch]
# names(Negative_set_seq_list_char_1000bp[[2]]) <- paste("neg", c(1:length(Negative_set_seq_list_1000bp[[2]])), sep = "_")

Positive_set_seq_list_char_1000bp[[2]] <- as.character(Positive_set_seq_list_1000bp[[2]])
aamatch <- match(Positive_set_seq_list_1000bp[[2]], Positive_set_seq_list_1000bp[[1]])
names(Positive_set_seq_list_char_1000bp[[2]]) <- names(Positive_set_seq_list_char_1000bp[[1]])[aamatch]

# names(Positive_set_seq_list_char_1000bp[[2]]) <- paste("pos", c(1:length(Positive_set_seq_list_1000bp[[2]])), sep = "_")

Negative_set_seq_list_char_1000bp[[4]] <- as.character(Negative_set_seq_list_1000bp[[4]])
aamatch <- match(Negative_set_seq_list_1000bp[[4]], Negative_set_seq_list_1000bp[[2]])
names(Negative_set_seq_list_char_1000bp[[4]]) <- names(Negative_set_seq_list_char_1000bp[[2]])[aamatch]
#names(Negative_set_seq_list_char_1000bp[[4]]) <- paste("neg", c(1:length(Negative_set_seq_list_1000bp[[4]])), sep = "_")

Positive_set_seq_list_char_1000bp[[4]] <- as.character(Positive_set_seq_list_1000bp[[4]])
aamatch <- match(Positive_set_seq_list_1000bp[[4]], Positive_set_seq_list_1000bp[[2]])
names(Positive_set_seq_list_char_1000bp[[4]]) <- names(Positive_set_seq_list_char_1000bp[[2]])[aamatch]
# names(Positive_set_seq_list_char_1000bp[[4]]) <- paste("pos", c(1:length(Positive_set_seq_list_1000bp[[4]])), sep = "_")


################################################################################################################
# motifScore for the third set
# running on hal
save(list = c("Negative_set_seq_list_1000bp",
              "TF.motifs.Expanded_new_pseudo_t"),
     file = "Seq_3_plus_new_motifs.RData")

cat(c("#!/bin/bash\n"), file = "copy_motifscore_3.sh", append = F)
for(i in 1:length(my_sample_3)){
  cat(paste0("rsync ~/remote/MotifScore_neg_", my_sample_3[i], ".RData ",
             "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/3/Neg_motifScore/", "\n"),
      file = "copy_motifscore_3.sh", append = T)
}


################################################################################################################
Positive_set_seq_list_chopped_parameters_1000bp <- list()
Positive_set_seq_list_chopped_1000bp <- list()
Positive_set_seq_list_chopped_parameters_1000bp[[1]] <- list(Enhancer_seq = "Positive_set_seq_list_char_1000bp[[1]]",
                                                             motifScore_output = "Positive_set_motifscore_1000bp[[1]]",
                                                             EnhancerGR = "Positive_set_gr_list_1000bp[[1]]",
                                                             motifList = "TF.motifs.Expanded_pseudo_t",
                                                             no_thresh = "T",
                                                             LLR_thresh = "0",
                                                             piece_length = "1000",
                                                             step_size = "200",
                                                             normalize_by_length = "F",
                                                             pairwise_closeby_site_freq = "T",
                                                             close_thr = "25",
                                                             LLR2PVal_lists = "TF.motifs.Expanded_LLR_to_pval_pseudo",
                                                             pval_threshold = "0.0001",
                                                             .diff_max_LLR = "F")

Positive_set_seq_list_chopped_1000bp[[1]] <- EnhancerChopper(Enhancer_seq = Positive_set_seq_list_char_1000bp[[1]],
                                                             motifScore_output = Positive_set_motifscore_1000bp[[1]],
                                                             EnhancerGR = Positive_set_gr_list_1000bp[[1]],
                                                             motifList = TF.motifs.Expanded_pseudo_t,
                                                             no_thresh = T,
                                                             LLR_thresh = 0,
                                                             piece_length = 1000,
                                                             step_size = 200,
                                                             normalize_by_length = F,
                                                             pairwise_closeby_site_freq = T,
                                                             close_thr = 25,
                                                             LLR2PVal_lists = TF.motifs.Expanded_LLR_to_pval_pseudo,
                                                             pval_threshold = 0.0001,
                                                             .diff_max_LLR = F)

Negative_set_seq_list_chopped_parameters_1000bp <- list()
Negative_set_seq_list_chopped_1000bp <- list()
Negative_set_seq_list_chopped_parameters_1000bp[[1]] <- list(Enhancer_seq = "Negative_set_seq_list_char_1000bp[[1]]",
                                                             motifScore_output = "Negative_set_motifscore_1000bp[[1]]",
                                                             EnhancerGR = "Negative_set_gr_list_1000bp[[1]]",
                                                             motifList = "TF.motifs.Expanded_pseudo_t",
                                                             no_thresh = "T",
                                                             LLR_thresh = "0",
                                                             piece_length = "1000",
                                                             step_size = "200",
                                                             normalize_by_length = "F",
                                                             pairwise_closeby_site_freq = "T",
                                                             close_thr = "25",
                                                             LLR2PVal_lists = "TF.motifs.Expanded_LLR_to_pval_pseudo",
                                                             pval_threshold = "0.0001",
                                                             .diff_max_LLR = "F")

Negative_set_seq_list_chopped_1000bp[[1]] <- EnhancerChopper(Enhancer_seq = Negative_set_seq_list_char_1000bp[[1]],
                                                             motifScore_output = Negative_set_motifscore_1000bp[[1]],
                                                             EnhancerGR = Negative_set_gr_list_1000bp[[1]],
                                                             motifList = TF.motifs.Expanded_pseudo_t,
                                                             no_thresh = T,
                                                             LLR_thresh = 0,
                                                             piece_length = 1000,
                                                             step_size = 200,
                                                             normalize_by_length = F,
                                                             pairwise_closeby_site_freq = T,
                                                             close_thr = 25,
                                                             LLR2PVal_lists = TF.motifs.Expanded_LLR_to_pval_pseudo,
                                                             pval_threshold = 0.0001,
                                                             .diff_max_LLR = F)

################################################################################################################
# for the third set choose 2500 negatives at random from the ones that motifscore has been calced:
aafiles <- list.files(path = "~/remote/", pattern ="MotifScore_neg_*")
aafiles <- unlist(lapply(strsplit(aafiles, split = "\\."), "[[", 1))
aafiles <- as.integer(unlist(lapply(strsplit(aafiles, split = "_"), "[[", 3)))
my_sample_3 <- sample(aafiles, size = 2500, replace = F)
Negative_set_seq_list_1000bp[[3]] <- Negative_set_seq_list_1000bp[[3]][my_sample_3]
# read the  negatives from hal output
Negative_set_motifscore_1000bp <- list()
Negative_set_motifscore_1000bp[[3]] <- list()
for(i in 1:length(Negative_set_seq_list_1000bp[[3]])){
  print(i)
  load(paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/3/Neg_motifScore/MotifScore_neg_", my_sample_3[i], ".RData"))
  Negative_set_motifscore_1000bp[[3]][[i]] <- my_motif_score
  rm("my_motif_score")
}

# get motif scores for positive set for RARA_half motif
aaPositive_set_motifscore_1000bp <- list()
aaPositive_set_motifscore_1000bp[[3]] <- list()
for(i in 1:length(Positive_set_seq_list_1000bp[[3]])){
  print(paste("positive seq nu ", i, " out of ", length(Positive_set_seq_list_1000bp[[3]])))
  aaPositive_set_motifscore_1000bp[[3]][[i]] <- MotifScore(seq = as.character(Positive_set_seq_list_1000bp[[3]])[i],
                                                  bg = c(0.25, 0.25, 0.25, 0.25),
                                                  motifList = TF.motifs.Expanded_new_pseudo_t[length(TF.motifs.Expanded_new_pseudo_t)])

}
Negative_set_seq_list_char_1000bp[[3]] <- as.character(Negative_set_seq_list_1000bp[[3]])
names(Negative_set_seq_list_char_1000bp[[3]]) <- paste("neg", c(1:length(Negative_set_seq_list_1000bp[[3]])), sep = "_")

Positive_set_seq_list_char_1000bp[[3]] <- as.character(Positive_set_seq_list_1000bp[[3]])
names(Positive_set_seq_list_char_1000bp[[3]]) <- paste("pos", c(1:length(Positive_set_seq_list_1000bp[[3]])), sep = "_")


aaNegative_set_motifscore_1000bp <- Negative_set_motifscore_1000bp
load("MotifScore_Results_full_1000bp.RData")
length(Positive_set_motifscore_1000bp[[1]][[1]]$Reverse)
Positive_set_motifscore_1000bp[[3]] <- Positive_set_motifscore_1000bp[[1]]
aatfmatch <- match(names(TF.motifs.Expanded_new_pseudo_t)[1:26], names(TF.motifs.Expanded_pseudo_t))
for(i in 1:length(Positive_set_motifscore_1000bp[[1]])){
  Positive_set_motifscore_1000bp[[3]][[i]]$Forward <- Positive_set_motifscore_1000bp[[1]][[i]]$Forward[aatfmatch]
  Positive_set_motifscore_1000bp[[3]][[i]]$Forward[[length(Positive_set_motifscore_1000bp[[3]][[i]]$Forward) + 1]] <- aaPositive_set_motifscore_1000bp[[3]][[i]]$Forward$RARA_half
  names(Positive_set_motifscore_1000bp[[3]][[i]]$Forward)[length(Positive_set_motifscore_1000bp[[3]][[i]]$Forward)] <- names(TF.motifs.Expanded_new_pseudo_t)[length(TF.motifs.Expanded_new_pseudo_t)]
  Positive_set_motifscore_1000bp[[3]][[i]]$Reverse <- Positive_set_motifscore_1000bp[[1]][[i]]$Reverse[aatfmatch]
  Positive_set_motifscore_1000bp[[3]][[i]]$Reverse[[length(Positive_set_motifscore_1000bp[[3]][[i]]$Reverse) + 1]] <- aaPositive_set_motifscore_1000bp[[3]][[i]]$Reverse$RARA_half
  names(Positive_set_motifscore_1000bp[[3]][[i]]$Reverse)[length(Positive_set_motifscore_1000bp[[3]][[i]]$Reverse)] <- names(TF.motifs.Expanded_new_pseudo_t)[length(TF.motifs.Expanded_new_pseudo_t)]
}
################################################################################################################

Negative_set_seq_list_chopped_1000bp[[3]] <- EnhancerChopper(Enhancer_seq = Negative_set_seq_list_char_1000bp[[3]],
                                                             motifScore_output = Negative_set_motifscore_1000bp[[3]],
                                                             EnhancerGR = Negative_set_gr_list_1000bp[[3]],
                                                             motifList = TF.motifs.Expanded_new_pseudo_t,
                                                             no_thresh = T,
                                                             LLR_thresh = 0,
                                                             piece_length = 1000,
                                                             step_size = 200,
                                                             normalize_by_length = F,
                                                             pairwise_closeby_site_freq = T,
                                                             close_thr = 25,
                                                             LLR2PVal_lists = TF.motifs.Expanded_LLR_to_pval_pseudo_2,
                                                             pval_threshold = 0.0001,
                                                             .diff_max_LLR = F)
Positive_set_seq_list_chopped_1000bp[[3]] <- EnhancerChopper(Enhancer_seq = Positive_set_seq_list_char_1000bp[[3]],
                                                             motifScore_output = Positive_set_motifscore_1000bp[[3]],
                                                             EnhancerGR = Positive_set_gr_list_1000bp[[3]],
                                                             motifList = TF.motifs.Expanded_new_pseudo_t,
                                                             no_thresh = T,
                                                             LLR_thresh = 0,
                                                             piece_length = 1000,
                                                             step_size = 200,
                                                             normalize_by_length = F,
                                                             pairwise_closeby_site_freq = T,
                                                             close_thr = 25,
                                                             LLR2PVal_lists = TF.motifs.Expanded_LLR_to_pval_pseudo_2,
                                                             pval_threshold = 0.0001,
                                                             .diff_max_LLR = F)

save(list = c("Positive_set_motifscore_1000bp", "Negative_set_motifscore_1000bp"), file = "MotifScore_Results_full_1000bp.RData")

################################################################################################################
## mapping LLR to p-value
## Writing Each TF in a different file
# repeating with the pseudo count motifs
TF.motifs.Expanded_pseudo_count <- lapply(TF.motifs.Expanded_pseudo, PWMtoCount, remove.degen=F, times=10000, pseudo=0)
for(i in 1:length(TF.motifs.Expanded_pseudo)){
  MotifWriter(motif.List = TF.motifs.Expanded_pseudo_count[i],  pseudo = 0.001,
              output.File.Name = paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Motifs_to_pval_pseudo/",names(TF.motifs.Expanded_pseudo_count)[i]))
}
aa <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Motifs_to_pval_pseudo_2", pattern = "*.Filter", full.names = T)
TF.motifs.Expanded_LLR_to_pval_pseudo <- list()
for(i in 1:length(aa)){
  TF.motifs.Expanded_LLR_to_pval_pseudo[[i]] <- read.table(aa[i], header = F)
}
names(TF.motifs.Expanded_LLR_to_pval_pseudo) <- c(names(TF.motifs.Expanded_pseudo)[1:5], "ESR1_4", names(TF.motifs.Expanded_pseudo)[6:33])

# get the LLR2Pval for the "RARA_half" motif
MotifWriter(motif.List = TF.motifs.Expanded_new_pseudo_t[27],  pseudo = 0.001,
            output.File.Name = paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Motifs_to_pval_pseudo_RARAhalf/",names(TF.motifs.Expanded_new_pseudo_t)[27]))
TF.motifs.Expanded_LLR_to_pval_pseudo_RARAHalf <- read.table("E_RNA/Motifs_to_pval_pseudo_RARAhalf/RARA_half.Filter", header = F)

TF.motifs.Expanded_LLR_to_pval_pseudo_2 <- TF.motifs.Expanded_LLR_to_pval_pseudo[match(names(TF.motifs.Expanded_new_pseudo_t)[1:26], names(TF.motifs.Expanded_pseudo_count))]
TF.motifs.Expanded_LLR_to_pval_pseudo_2[[27]] <- TF.motifs.Expanded_LLR_to_pval_pseudo_RARAHalf
names(TF.motifs.Expanded_LLR_to_pval_pseudo_2)[27] <- names(TF.motifs.Expanded_new_pseudo_t)[27]
################################################################################################################
# Set some maximum length and chop enhancers above that length
# chopping with finding pairwise site frequency
Positive_set_seq_list_char <- list()
Positive_set_seq_list_char[[1]] <- as.character(Positive_set_seq_list[[1]])
Negative_set_seq_list_char <- list()
Negative_set_seq_list_char[[1]] <- as.character(Negative_set_seq_list[[1]])
names(Positive_set_seq_list_char[[1]]) <- paste0("pos_", c(1:length(Positive_set_seq_list_char[[1]])))
names(Negative_set_seq_list_char[[1]]) <- paste0("neg_", c(1:length(Negative_set_seq_list_char[[1]])))

Positive_set_seq_list_chopped <- list()
Positive_set_seq_list_chopped_parameters <- list()
Positive_set_seq_list_chopped_parameters[[1]] <- list(Enhancer_seq = "Positive_set_seq_list_char[[1]]",
                                                      motifScore_output = "Positive_set_motifscore[[1]]",
                                                      EnhancerGR = "Positive_set_gr_list[[1]]",
                                                      motifList = "TF.motifs.Expanded_pseudo_t",
                                                      no_thresh = "T",
                                                      LLR_thresh = "0",
                                                      piece_length = "1350",
                                                      step_size = "200",
                                                      normalize_by_length = "F",
                                                      pairwise_closeby_site_freq = "T",
                                                      close_thr = "25",
                                                      LLR2PVal_lists = "TF.motifs.Expanded_LLR_to_pval_pseudo",
                                                      pval_threshold = "0.0001",
                                                      .diff_max_LLR = "F")

Positive_set_seq_list_chopped[[1]] <- EnhancerChopper(Enhancer_seq = Positive_set_seq_list_char[[1]],
                                                      motifScore_output = Positive_set_motifscore[[1]],
                                                      EnhancerGR = Positive_set_gr_list[[1]],
                                                      motifList = TF.motifs.Expanded_pseudo_t,
                                                      no_thresh = T,
                                                      LLR_thresh = 0,
                                                      piece_length = 1350,
                                                      step_size = 200,
                                                      normalize_by_length = F,
                                                      pairwise_closeby_site_freq = T,
                                                      close_thr = 25,
                                                      LLR2PVal_lists = TF.motifs.Expanded_LLR_to_pval_pseudo,
                                                      pval_threshold = 0.0001,
                                                      .diff_max_LLR = F)

Negative_set_seq_list_chopped <- list()
Negative_set_seq_list_chopped_parameters <- list()
Negative_set_seq_list_chopped_parameters[[1]] <- list(Enhancer_seq = "Negative_set_seq_list_char[[1]]",
                                                      motifScore_output = "Negative_set_motifscore[[1]]",
                                                      EnhancerGR = "Negative_set_gr_list[[1]]",
                                                      motifList = "TF.motifs.Expanded_pseudo_t",
                                                      no_thresh = "T",
                                                      LLR_thresh = "0",
                                                      piece_length = "1350",
                                                      step_size = "200",
                                                      normalize_by_length = "F",
                                                      pairwise_closeby_site_freq = "T",
                                                      close_thr = "25",
                                                      LLR2PVal_lists = "TF.motifs.Expanded_LLR_to_pval_pseudo",
                                                      pval_threshold = "0.0001",
                                                      .diff_max_LLR = "F")
Negative_set_seq_list_chopped[[1]] <- EnhancerChopper(Enhancer_seq = Negative_set_seq_list_char[[1]],
                                                      motifScore_output = Negative_set_motifscore[[1]],
                                                      EnhancerGR = Negative_set_gr_list[[1]],
                                                      motifList = TF.motifs.Expanded_pseudo_t,
                                                      no_thresh = T,
                                                      LLR_thresh = 0,
                                                      piece_length = 1350,
                                                      step_size = 200,
                                                      normalize_by_length = F,
                                                      pairwise_closeby_site_freq = T,
                                                      close_thr = 25,
                                                      LLR2PVal_lists = TF.motifs.Expanded_LLR_to_pval_pseudo,
                                                      pval_threshold = 0.0001,
                                                      .diff_max_LLR = F)
# concatenate all scores in a matrix where rows are name by sequences also noting the length of each sequence
Negative_set_seq_list_chopped_scoreMat_list <- list() 
Negative_set_seq_list_chopped_adjanScore_list <- list() 
Negative_set_seq_list_chopped_OverlapScore_list <- list() 

Negative_set_seq_list_chopped_scoreMat_list[[1]] <- do.call(rbind, Negative_set_seq_list_chopped[[1]]$Chopped_Score)
Negative_set_seq_list_chopped_scoreMat_list[[1]] <- cbind(Negative_set_seq_list_chopped_scoreMat_list[[1]], 
                                                          nchar(unlist(Negative_set_seq_list_chopped[[1]]$Chopped_Seq)))
colnames(Negative_set_seq_list_chopped_scoreMat_list[[1]])[ncol(Negative_set_seq_list_chopped_scoreMat_list[[1]])] <- "length"
head(Negative_set_seq_list_chopped_scoreMat_list[[1]])

Negative_set_seq_list_chopped_adjanScore_list[[1]] <- do.call(rbind, Negative_set_seq_list_chopped[[1]]$Adjan_list)
Negative_set_seq_list_chopped_OverlapScore_list[[1]] <- do.call(rbind, Negative_set_seq_list_chopped[[1]]$overlap_list)

Positive_set_seq_list_chopped_scoreMat_list <- list() 
Positive_set_seq_list_chopped_adjanScore_list <- list() 
Positive_set_seq_list_chopped_OverlapScore_list <- list() 

Positive_set_seq_list_chopped_scoreMat_list[[1]] <- do.call(rbind, Positive_set_seq_list_chopped[[1]]$Chopped_Score)
Positive_set_seq_list_chopped_scoreMat_list[[1]] <- cbind(Positive_set_seq_list_chopped_scoreMat_list[[1]], 
                                                          nchar(unlist(Positive_set_seq_list_chopped[[1]]$Chopped_Seq)))
colnames(Positive_set_seq_list_chopped_scoreMat_list[[1]])[ncol(Positive_set_seq_list_chopped_scoreMat_list[[1]])] <- "length"
head(Positive_set_seq_list_chopped_scoreMat_list[[1]])

Positive_set_seq_list_chopped_adjanScore_list[[1]] <- do.call(rbind, Positive_set_seq_list_chopped[[1]]$Adjan_list)
Positive_set_seq_list_chopped_OverlapScore_list[[1]] <- do.call(rbind, Positive_set_seq_list_chopped[[1]]$overlap_list)



par(mfrow = c(2, 1), mar = c(4,4,4,4))
boxplot.matrix(Positive_set_seq_list_chopped_scoreMat_list[[1]][,1:33], ylim = c(0, 65000),
               las = 2, main = "positive_Set_sum_LR", outline = F)
abline(h = seq(0, 60000, 5000), col = 2, lty = 4, lwd = 0.5)
boxplot.matrix(Negative_set_seq_list_chopped_scoreMat_list[[1]][,1:33],ylim = c(0, 65000),
               las = 2, main = "negative_Set_sum_LR", outline = F)
abline(h = seq(0, 60000, 5000), col = 2, lty = 4, lwd = 0.5)

aa_newmat <- matrix(nrow = max(nrow(Positive_set_seq_list_chopped_scoreMat_list[[1]]),
                               nrow(Negative_set_seq_list_chopped_scoreMat_list[[1]])),
                    ncol = 66)
colnames(aa_newmat) <- character(66)
for(i in 1:(ncol(Positive_set_seq_list_chopped_scoreMat_list[[1]]) - 1)){
  aa_newmat[1:nrow(Positive_set_seq_list_chopped_scoreMat_list[[1]]),(2*i - 1)] <- Positive_set_seq_list_chopped_scoreMat_list[[1]][, i]
  colnames(aa_newmat)[2*i - 1] <- paste(colnames(Positive_set_seq_list_chopped_scoreMat_list[[1]])[i], "Pos", sep = "_")
  aa_newmat[1:nrow(Negative_set_seq_list_chopped_scoreMat_list[[1]]),(2*i)] <- Negative_set_seq_list_chopped_scoreMat_list[[1]][, i]
  colnames(aa_newmat)[2*i] <- paste(colnames(Negative_set_seq_list_chopped_scoreMat_list[[1]])[i], "Neg", sep = "_")
  
}
par(mfrow = c(1, 1), mar = c(6,4,4,4))

boxplot.matrix(aa_newmat,
               las = 2, main = "sum_LLR", outline = F, col = rep(c("green", "red"), 33))
legend("top", legend = c("positive", "negative"), fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)

par(mfrow = c(1, 1), mar = c(6,4,4,4))

boxplot.matrix(aa_newmat,
               las = 2, main = "sum_LLR", outline = F, col = rep(c("green", "red"), 33), ylim = c(0, 20000))
legend("top", legend = c("positive", "negative"), fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)

###
# with normalizing for length
Positive_set_seq_list_chopped_parameters[[2]] <- list(Enhancer_seq = "Positive_set_seq_list_char[[1]]",
                                                      motifScore_output = "Positive_set_motifscore[[1]]",
                                                      EnhancerGR = "Positive_set_gr_list[[1]]",
                                                      motifList = "TF.motifs.Expanded_pseudo_t",
                                                      no_thresh = "T",
                                                      LLR_thresh = "0",
                                                      piece_length = "1350",
                                                      step_size = "200",
                                                      normalize_by_length = "T",
                                                      pairwise_closeby_site_freq = "T",
                                                      close_thr = "25",
                                                      LLR2PVal_lists = "TF.motifs.Expanded_LLR_to_pval_pseudo",
                                                      pval_threshold = "0.0001",
                                                      .diff_max_LLR = "F")
  
Positive_set_seq_list_chopped[[2]] <- EnhancerChopper(Enhancer_seq = Positive_set_seq_list_char[[1]],
                                                      motifScore_output = Positive_set_motifscore[[1]],
                                                      EnhancerGR = Positive_set_gr_list[[1]],
                                                      motifList = TF.motifs.Expanded_pseudo_t,
                                                      no_thresh = T,
                                                      LLR_thresh = 0,
                                                      piece_length = 1350,
                                                      step_size = 200,
                                                      normalize_by_length = T,
                                                      pairwise_closeby_site_freq = T,
                                                      close_thr = 25,
                                                      LLR2PVal_lists = TF.motifs.Expanded_LLR_to_pval_pseudo,
                                                      pval_threshold = 0.0001,
                                                      .diff_max_LLR = F)
Negative_set_seq_list_chopped_parameters[[2]] <- list(Enhancer_seq = "Negative_set_seq_list_char[[1]]",
                                                      motifScore_output = "Negative_set_motifscore[[1]]",
                                                      EnhancerGR = "Negative_set_gr_list[[1]]",
                                                      motifList = "TF.motifs.Expanded_pseudo_t",
                                                      no_thresh = "T",
                                                      LLR_thresh = "0",
                                                      piece_length = "1350",
                                                      step_size = "200",
                                                      normalize_by_length = "T",
                                                      pairwise_closeby_site_freq = "T",
                                                      close_thr = "25",
                                                      LLR2PVal_lists = "TF.motifs.Expanded_LLR_to_pval_pseudo",
                                                      pval_threshold = "0.0001",
                                                      .diff_max_LLR = "F")

Negative_set_seq_list_chopped[[2]] <- EnhancerChopper(Enhancer_seq = Negative_set_seq_list_char[[1]],
                                                      motifScore_output = Negative_set_motifscore[[1]],
                                                      EnhancerGR = Negative_set_gr_list[[1]],
                                                      motifList = TF.motifs.Expanded_pseudo_t,
                                                      no_thresh = T,
                                                      LLR_thresh = 0,
                                                      piece_length = 1350,
                                                      step_size = 200,
                                                      normalize_by_length = T, 
                                                      pairwise_closeby_site_freq = T,
                                                      close_thr = 25,
                                                      LLR2PVal_lists = TF.motifs.Expanded_LLR_to_pval_pseudo,
                                                      pval_threshold = 0.0001,
                                                      .diff_max_LLR = F)


Negative_set_seq_list_chopped_scoreMat_list[[2]] <- do.call(rbind, Negative_set_seq_list_chopped[[2]]$Chopped_Score)
Negative_set_seq_list_chopped_scoreMat_list[[2]] <- cbind(Negative_set_seq_list_chopped_scoreMat_list[[2]], 
                                                          nchar(unlist(Negative_set_seq_list_chopped[[2]]$Chopped_Seq)))
colnames(Negative_set_seq_list_chopped_scoreMat_list[[2]])[ncol(Negative_set_seq_list_chopped_scoreMat_list[[2]])] <- "length"
head(Negative_set_seq_list_chopped_scoreMat_list[[2]])

Negative_set_seq_list_chopped_adjanScore_list[[2]] <- do.call(rbind, Negative_set_seq_list_chopped[[2]]$Adjan_list)
Negative_set_seq_list_chopped_OverlapScore_list[[2]] <- do.call(rbind, Negative_set_seq_list_chopped[[2]]$overlap_list)


Positive_set_seq_list_chopped_scoreMat_list[[2]] <- do.call(rbind, Positive_set_seq_list_chopped[[2]]$Chopped_Score)
Positive_set_seq_list_chopped_scoreMat_list[[2]] <- cbind(Positive_set_seq_list_chopped_scoreMat_list[[2]], 
                                                          nchar(unlist(Positive_set_seq_list_chopped[[2]]$Chopped_Seq)))
colnames(Positive_set_seq_list_chopped_scoreMat_list[[2]])[ncol(Positive_set_seq_list_chopped_scoreMat_list[[2]])] <- "length"
head(Positive_set_seq_list_chopped_scoreMat_list[[2]])

Positive_set_seq_list_chopped_adjanScore_list[[2]] <- do.call(rbind, Positive_set_seq_list_chopped[[2]]$Adjan_list)
Positive_set_seq_list_chopped_OverlapScore_list[[2]] <- do.call(rbind, Positive_set_seq_list_chopped[[2]]$overlap_list)


par(mfrow = c(2, 1), mar = c(4,4,4,4))
boxplot.matrix(Positive_set_seq_list_chopped_scoreMat_list[[2]][,1:33], ylim = c(0, 150000),
               las = 2, main = "positive_Set_sum_LR", outline = F)
abline(h = seq(0, 150000, 5000), col = 2, lty = 4, lwd = 0.5)
boxplot.matrix(Negative_set_seq_list_chopped_scoreMat_list[[2]][,1:33],ylim = c(0, 150000),
               las = 2, main = "negative_Set_sum_LR", outline = F)
abline(h = seq(0, 150000, 5000), col = 2, lty = 4, lwd = 0.5)

aa_newmat <- matrix(nrow = max(nrow(Positive_set_seq_list_chopped_scoreMat_list[[2]]),
                               nrow(Negative_set_seq_list_chopped_scoreMat_list[[2]])),
                    ncol = 66)
colnames(aa_newmat) <- character(66)
for(i in 1:(ncol(Positive_set_seq_list_chopped_scoreMat_list[[2]]) - 1)){
  aa_newmat[1:nrow(Positive_set_seq_list_chopped_scoreMat_list[[2]]),(2*i - 1)] <- Positive_set_seq_list_chopped_scoreMat_list[[2]][, i]
  colnames(aa_newmat)[2*i - 1] <- paste(colnames(Positive_set_seq_list_chopped_scoreMat_list[[2]])[i], "Pos", sep = "_")
  aa_newmat[1:nrow(Negative_set_seq_list_chopped_scoreMat_list[[2]]),(2*i)] <- Negative_set_seq_list_chopped_scoreMat_list[[2]][, i]
  colnames(aa_newmat)[2*i] <- paste(colnames(Negative_set_seq_list_chopped_scoreMat_list[[2]])[i], "Neg", sep = "_")
  
}

par(mfrow = c(1, 1), mar = c(6,4,4,4))

boxplot.matrix(aa_newmat,
               las = 2, main = "sum_LLR_length_normalized", outline = F, col = rep(c("green", "red"), 33))
legend("top", legend = c("positive", "negative"), fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)

par(mfrow = c(1, 1), mar = c(6,4,4,4))

boxplot.matrix(aa_newmat,
               las = 2, main = "sum_LLR_length_normalized", outline = F, col = rep(c("green", "red"), 33), ylim = c(0, 20000))
legend("top", legend = c("positive", "negative"), fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)


aa_newmat2 <- matrix(nrow = max(nrow(Positive_set_seq_list_chopped_adjanScore_list[[2]]),
                               nrow(Negative_set_seq_list_chopped_adjanScore_list[[2]])),
                    ncol = sum(ncol(Positive_set_seq_list_chopped_adjanScore_list[[2]]),
                               ncol(Negative_set_seq_list_chopped_adjanScore_list[[2]])))
colnames(aa_newmat2) <- character(sum(ncol(Positive_set_seq_list_chopped_adjanScore_list[[2]]),
                                      ncol(Positive_set_seq_list_chopped_adjanScore_list[[2]])))
for(i in 1:ncol(Positive_set_seq_list_chopped_adjanScore_list[[2]])){
  aa_newmat2[1:nrow(Positive_set_seq_list_chopped_adjanScore_list[[2]]),(2*i - 1)] <- Positive_set_seq_list_chopped_adjanScore_list[[2]][, i]
  colnames(aa_newmat2)[2*i - 1] <- paste(colnames(Positive_set_seq_list_chopped_adjanScore_list[[2]])[i], "Pos", sep = "_")
  aa_newmat2[1:nrow(Negative_set_seq_list_chopped_adjanScore_list[[2]]),(2*i)] <- Negative_set_seq_list_chopped_adjanScore_list[[2]][, i]
  colnames(aa_newmat2)[2*i] <- paste(colnames(Negative_set_seq_list_chopped_adjanScore_list[[2]])[i], "Neg", sep = "_")
  
}
par(mfrow = c(1, 1), mar = c(15,4,4,4))

boxplot.matrix(aa_newmat2,
               las = 2, main = "adj score", outline = F, col = rep(c("green", "red"), ncol(Positive_set_seq_list_chopped_adjanScore_list[[2]])))
legend("top", legend = c("positive", "negative"), fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)

jpeg(filename = "adj_scores_pos_neg.jpeg",
     width = 400, height = 10, units = "cm", res = 600)
boxplot.matrix(aa_newmat2,
               las = 2, main = "adj score", outline = F, col = rep(c("green", "red"), ncol(Positive_set_seq_list_chopped_adjanScore_list[[2]])))
dev.off()

#############################################################################
# changing the motif set used for MotifScanning
aan <- which(names(TF.motifs.Expanded_pseudo_t) %in% c("NFIB", "NR3C1", "PAX2", "RARA", "RARG", "RXRA", "RXRB"))

i <- 6
names(TF.motifs.Expanded_pseudo_t[aan[i]])
seqLogo::seqLogo(TF.motifs.Expanded_pseudo_t[[33]])
TF.motifs.Expanded_new_pseudo_t <- TF.motifs.Expanded_pseudo_t[-aan]
TF.motifs.Expanded_new_pseudo_t[[length(TF.motifs.Expanded_new_pseudo_t) + 1]]  <- TF.motifs.Expanded_pseudo_t$RARA[,1:6]
names(TF.motifs.Expanded_new_pseudo_t)[length(TF.motifs.Expanded_new_pseudo_t)] <- "RARA_half"

#############################################################################
# concatenate all scores in a matrix where rows are named by sequences --> for the 1000 bp version
Negative_set_seq_list_chopped_scoreMat_list_1000bp <- list() 
Negative_set_seq_list_chopped_adjanScore_list_1000bp <- list() 
Negative_set_seq_list_chopped_OverlapScore_list_1000bp <- list() 

Negative_set_seq_list_chopped_scoreMat_list_1000bp[[1]] <- do.call(rbind, Negative_set_seq_list_chopped_1000bp[[1]]$Chopped_Score)
head(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[1]])

Negative_set_seq_list_chopped_adjanScore_list_1000bp[[1]] <- do.call(rbind, Negative_set_seq_list_chopped_1000bp[[1]]$Adjan_list)
Negative_set_seq_list_chopped_OverlapScore_list_1000bp[[1]] <- do.call(rbind, Negative_set_seq_list_chopped_1000bp[[1]]$overlap_list)

Positive_set_seq_list_chopped_scoreMat_list_1000bp <- list() 
Positive_set_seq_list_chopped_adjanScore_list_1000bp <- list() 
Positive_set_seq_list_chopped_OverlapScore_list_1000bp <- list() 

Positive_set_seq_list_chopped_scoreMat_list_1000bp[[1]] <- do.call(rbind, Positive_set_seq_list_chopped_1000bp[[1]]$Chopped_Score)
head(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[1]])

Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]] <- do.call(rbind, Positive_set_seq_list_chopped_1000bp[[1]]$Adjan_list)
Positive_set_seq_list_chopped_OverlapScore_list_1000bp[[1]] <- do.call(rbind, Positive_set_seq_list_chopped_1000bp[[1]]$overlap_list)


par(mfrow = c(2, 1), mar = c(4,4,4,4))
boxplot.matrix(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[1]],
               ylim = c(0, 30000),
               las = 2, main = "positive_Set_sum_LR_1000bp", outline = F)
abline(h = seq(0, 150000, 5000), col = 2, lty = 4, lwd = 0.5)
boxplot.matrix(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[1]],
               ylim = c(0, 30000),
               las = 2, main = "negative_Set_sum_LR_1000bp", outline = F)
abline(h = seq(0, 150000, 5000), col = 2, lty = 4, lwd = 0.5)

aa_newmat <- matrix(nrow = max(nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[1]]),
                               nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[1]])),
                    ncol = 66)
colnames(aa_newmat) <- character(66)
for(i in 1:ncol(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[1]])){
  aa_newmat[1:nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[1]]),(2*i - 1)] <- Positive_set_seq_list_chopped_scoreMat_list_1000bp[[1]][, i]
  colnames(aa_newmat)[2*i - 1] <- paste(colnames(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[1]])[i], "Pos", sep = "_")
  aa_newmat[1:nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[1]]),(2*i)] <- Negative_set_seq_list_chopped_scoreMat_list_1000bp[[1]][, i]
  colnames(aa_newmat)[2*i] <- paste(colnames(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[1]])[i], "Neg", sep = "_")
  
}

par(mfrow = c(1, 1), mar = c(6,4,4,4))

boxplot.matrix(aa_newmat,
               las = 2, main = "sum_LLR_1000bp", outline = F, col = rep(c("green", "red"), 33), ylim = c(0, 2000))
legend("top", legend = c("positive", "negative"), fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)



aa_newmat2 <- matrix(nrow = max(nrow(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]]),
                                nrow(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[1]])),
                     ncol = sum(ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]]),
                                ncol(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[1]])))
colnames(aa_newmat2) <- character(sum(ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]]),
                                      ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]])))
for(i in 1:ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]])){
  aa_newmat2[1:nrow(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]]),(2*i - 1)] <- Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]][, i]
  colnames(aa_newmat2)[2*i - 1] <- paste(colnames(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]])[i], "Pos", sep = "_")
  aa_newmat2[1:nrow(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[1]]),(2*i)] <- Negative_set_seq_list_chopped_adjanScore_list_1000bp[[1]][, i]
  colnames(aa_newmat2)[2*i] <- paste(colnames(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[1]])[i], "Neg", sep = "_")
  
}
par(mfrow = c(1, 1), mar = c(15,4,4,4))

boxplot.matrix(aa_newmat2,
               las = 2, main = "adj score", outline = F, col = rep(c("green", "red"), ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]])))
legend("top", legend = c("positive", "negative"), fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)

jpeg(filename = "adj_scores_pos_neg_1000bp.jpeg",
     width = 400, height = 10, units = "cm", res = 600)
boxplot.matrix(aa_newmat2,
               las = 2, main = "adj score", outline = F,
               col = rep(c("green", "red"),
                         ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]])))
dev.off()
#############################################################################
# creating the above matrices for the second positive_negative set of 1000bp length
aamatch <- match(Negative_set_seq_list_1000bp[[2]], Negative_set_seq_list_1000bp[[1]])
all(Negative_set_seq_list_1000bp[[1]][aamatch] == Negative_set_seq_list_1000bp[[2]])


Negative_set_seq_list_chopped_scoreMat_list_1000bp[[2]] <- Negative_set_seq_list_chopped_scoreMat_list_1000bp[[1]][aamatch,]
head(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[2]])
Negative_set_seq_list_chopped_adjanScore_list_1000bp[[2]] <- Negative_set_seq_list_chopped_adjanScore_list_1000bp[[1]][aamatch,]
Negative_set_seq_list_chopped_OverlapScore_list_1000bp[[2]] <- Negative_set_seq_list_chopped_OverlapScore_list_1000bp[[1]][aamatch,]

Positive_set_seq_list_chopped_scoreMat_list_1000bp[[2]] <- Positive_set_seq_list_chopped_scoreMat_list_1000bp[[1]]
Positive_set_seq_list_chopped_adjanScore_list_1000bp[[2]] <- Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]]
Positive_set_seq_list_chopped_OverlapScore_list_1000bp[[2]] <- Positive_set_seq_list_chopped_OverlapScore_list_1000bp[[1]]


par(mfrow = c(2, 1), mar = c(4,4,4,4))
boxplot.matrix(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[2]],
               ylim = c(0, 30000),
               las = 2, main = "positive_Set_sum_LR_1000bp", outline = F)
abline(h = seq(0, 150000, 5000), col = 2, lty = 4, lwd = 0.5)
boxplot.matrix(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[2]],
               ylim = c(0, 30000),
               las = 2, main = "negative_Set_sum_LR_1000bp_v2", outline = F)
abline(h = seq(0, 150000, 5000), col = 2, lty = 4, lwd = 0.5)

aa_newmat <- matrix(nrow = max(nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[2]]),
                               nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[2]])),
                    ncol = 66)
colnames(aa_newmat) <- character(66)
for(i in 1:ncol(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[2]])){
  aa_newmat[1:nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[2]]),(2*i - 1)] <- Positive_set_seq_list_chopped_scoreMat_list_1000bp[[2]][, i]
  colnames(aa_newmat)[2*i - 1] <- paste(colnames(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[2]])[i], "Pos", sep = "_")
  aa_newmat[1:nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[2]]),(2*i)] <- Negative_set_seq_list_chopped_scoreMat_list_1000bp[[2]][, i]
  colnames(aa_newmat)[2*i] <- paste(colnames(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[2]])[i], "Neg", sep = "_")
  
}

par(mfrow = c(1, 1), mar = c(6,4,4,4))

boxplot.matrix(aa_newmat,
               las = 2, main = "sum_LLR_1000bp", outline = F, col = rep(c("green", "red"), 33), ylim = c(0,5000))
legend("top", legend = c("positive", "negative"), fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)

summary(aa_newmat["NR3C1_Pos"], is.)

aa_newmat2 <- matrix(nrow = max(nrow(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[2]]),
                                nrow(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[2]])),
                     ncol = sum(ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[2]]),
                                ncol(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[2]])))
colnames(aa_newmat2) <- character(sum(ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[2]]),
                                      ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[2]])))
for(i in 1:ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[2]])){
  aa_newmat2[1:nrow(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[2]]),(2*i - 1)] <- Positive_set_seq_list_chopped_adjanScore_list_1000bp[[2]][, i]
  colnames(aa_newmat2)[2*i - 1] <- paste(colnames(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[2]])[i], "Pos", sep = "_")
  aa_newmat2[1:nrow(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[2]]),(2*i)] <- Negative_set_seq_list_chopped_adjanScore_list_1000bp[[2]][, i]
  colnames(aa_newmat2)[2*i] <- paste(colnames(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[2]])[i], "Neg", sep = "_")
  
}
par(mfrow = c(1, 1), mar = c(15,4,4,4))

boxplot.matrix(aa_newmat2,
               las = 2, main = "adj score", outline = F, col = rep(c("green", "red"), ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]])))
legend("top", legend = c("positive", "negative"), fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)

jpeg(filename = "adj_scores_pos_neg_1000bp.jpeg",
     width = 400, height = 10, units = "cm", res = 600)
boxplot.matrix(aa_newmat2,
               las = 2, main = "adj score", outline = F,
               col = rep(c("green", "red"),
                         ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]])))
dev.off()
#############################################################################
# for the third neg_pos set
Negative_set_seq_list_chopped_scoreMat_list_1000bp[[3]] <- do.call(rbind, Negative_set_seq_list_chopped_1000bp[[3]]$Chopped_Score)
head(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[3]])

Negative_set_seq_list_chopped_adjanScore_list_1000bp[[3]] <- do.call(rbind, Negative_set_seq_list_chopped_1000bp[[3]]$Adjan_list)
Negative_set_seq_list_chopped_OverlapScore_list_1000bp[[3]] <- do.call(rbind, Negative_set_seq_list_chopped_1000bp[[3]]$overlap_list)


Positive_set_seq_list_chopped_scoreMat_list_1000bp[[3]] <-do.call(rbind, Positive_set_seq_list_chopped_1000bp[[3]]$Chopped_Score)
Positive_set_seq_list_chopped_adjanScore_list_1000bp[[3]] <- do.call(rbind, Positive_set_seq_list_chopped_1000bp[[3]]$Adjan_list)
Positive_set_seq_list_chopped_OverlapScore_list_1000bp[[3]] <- do.call(rbind, Positive_set_seq_list_chopped_1000bp[[3]]$overlap_list)


par(mfrow = c(2, 1), mar = c(4,4,4,4))
boxplot.matrix(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[3]],
               ylim = c(0, 30000),
               las = 2, main = "positive_Set_sum_LR_1000bp", outline = F)
abline(h = seq(0, 150000, 5000), col = 2, lty = 4, lwd = 0.5)
boxplot.matrix(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[3]],
               ylim = c(0, 30000),
               las = 2, main = "negative_Set_sum_LR_1000bp", outline = F)
abline(h = seq(0, 150000, 5000), col = 2, lty = 4, lwd = 0.5)

aa_newmat <- matrix(nrow = max(nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[3]]),
                               nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[3]])),
                    ncol = 2 * ncol(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[3]]))
colnames(aa_newmat) <- character(ncol(aa_newmat))
for(i in 1:ncol(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[3]])){
  aa_newmat[1:nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[3]]),(2*i - 1)] <- Positive_set_seq_list_chopped_scoreMat_list_1000bp[[3]][, i]
  colnames(aa_newmat)[2*i - 1] <- paste(colnames(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[3]])[i], "Pos", sep = "_")
  aa_newmat[1:nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[3]]),(2*i)] <- Negative_set_seq_list_chopped_scoreMat_list_1000bp[[3]][, i]
  colnames(aa_newmat)[2*i] <- paste(colnames(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[3]])[i], "Neg", sep = "_")
  
}

par(mfrow = c(1, 1), mar = c(6,4,4,4))

boxplot.matrix(aa_newmat,
               las = 2, main = "sum_LLR_1000bp", outline = F, col = rep(c("green", "red"), 33), ylim = c(0, 20000))
legend("top", legend = c("positive", "negative"), fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)
#############################################################################
#############################################################################
# creating same matrices for set 4
aamatch <- match(Negative_set_seq_list_1000bp[[4]], Negative_set_seq_list_1000bp[[2]])
all(Negative_set_seq_list_1000bp[[2]][aamatch] == Negative_set_seq_list_1000bp[[4]])


Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]] <- Negative_set_seq_list_chopped_scoreMat_list_1000bp[[2]][aamatch,]
head(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]])
Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]] <- Negative_set_seq_list_chopped_adjanScore_list_1000bp[[2]][aamatch,]
Negative_set_seq_list_chopped_OverlapScore_list_1000bp[[4]] <- Negative_set_seq_list_chopped_OverlapScore_list_1000bp[[2]][aamatch,]

aamatch <- match(Positive_set_seq_list_1000bp[[4]], Positive_set_seq_list_1000bp[[2]])
all(Positive_set_seq_list_1000bp[[2]][aamatch] == Positive_set_seq_list_1000bp[[4]])

Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]] <- Positive_set_seq_list_chopped_scoreMat_list_1000bp[[2]][aamatch,]
Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]] <- Positive_set_seq_list_chopped_adjanScore_list_1000bp[[2]][aamatch,]
Positive_set_seq_list_chopped_OverlapScore_list_1000bp[[4]] <- Positive_set_seq_list_chopped_OverlapScore_list_1000bp[[2]][aamatch,]


par(mfrow = c(2, 1), mar = c(4,4,4,4))
boxplot.matrix(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]],
               ylim = c(0, 30000),
               las = 2, main = "positive_Set4_sum_LR_1000bp", outline = F)
abline(h = seq(0, 150000, 5000), col = 2, lty = 4, lwd = 0.5)
boxplot.matrix(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]],
               ylim = c(0, 30000),
               las = 2, main = "negative_Set4_sum_LR_1000bp", outline = F)
abline(h = seq(0, 150000, 5000), col = 2, lty = 4, lwd = 0.5)

aa_newmat <- matrix(nrow = max(nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]]),
                               nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]])),
                    ncol = 66)
colnames(aa_newmat) <- character(66)
for(i in 1:ncol(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]])){
  aa_newmat[1:nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]]),(2*i - 1)] <- Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]][, i]
  colnames(aa_newmat)[2*i - 1] <- paste(colnames(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]])[i], "Pos", sep = "_")
  aa_newmat[1:nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]]),(2*i)] <- Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]][, i]
  colnames(aa_newmat)[2*i] <- paste(colnames(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]])[i], "Neg", sep = "_")
  
}

par(mfrow = c(1, 1), mar = c(6,4,4,4))

boxplot.matrix(aa_newmat,
               las = 2, main = "sum_LLR_1000bp_4", outline = F, col = rep(c("green", "red"), 33), ylim = c(0,20000))
legend("top", legend = c("positive", "negative"), fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)


aa_newmat2 <- matrix(nrow = max(nrow(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]]),
                                nrow(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]])),
                     ncol = sum(ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]]),
                                ncol(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]])))
colnames(aa_newmat2) <- character(sum(ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]]),
                                      ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]])))
for(i in 1:ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]])){
  aa_newmat2[1:nrow(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]]),(2*i - 1)] <- Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]][, i]
  colnames(aa_newmat2)[2*i - 1] <- paste(colnames(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]])[i], "Pos", sep = "_")
  aa_newmat2[1:nrow(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]]),(2*i)] <- Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]][, i]
  colnames(aa_newmat2)[2*i] <- paste(colnames(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]])[i], "Neg", sep = "_")
  
}
par(mfrow = c(1, 1), mar = c(15,4,4,4))

boxplot.matrix(aa_newmat2,
               las = 2, main = "adj score", outline = F, col = rep(c("green", "red"), ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]])))
legend("top", legend = c("positive", "negative"), fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)

jpeg(filename = "adj_scores_pos_neg_1000bp_set4.jpeg",
     width = 400, height = 10, units = "cm", res = 600)
boxplot.matrix(aa_newmat2,
               las = 2, main = "adj score", outline = F,
               col = rep(c("green", "red"),
                         ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]])))
dev.off()



#############################################################################
# writing jobs to be submitted to hal
#cat(c("#!/bin/bash\n"), file = "run_motifscore.run", append = F)
for(i in 570:2349){
  cat(paste0("Rscript MotifScoreCalc_neg.R ",
             i,
             "\n"),
      file = "run_motifscore.run", append = T)
}

for(i in 1:length(Negative_set_seq_list_1000bp[[3]])){
  cat(paste0("Rscript MotifScoreCalc_neg.R ",
             i,
             "\n"),
      file = "run_motifscore_3.run", append = T)
}


#############################################################################
################################################################################################################


GEMSTAT_Ensemble_train_SeqList[[2]]

aamatch <- match(Negative_set_seq_list_1000bp[[2]], Negative_set_seq_list_1000bp[[1]])
all(Negative_set_seq_list_1000bp[[1]][aamatch] == Negative_set_seq_list_1000bp[[2]])


Negative_set_seq_list_chopped_scoreMat_list_1000bp[[2]] <- Negative_set_seq_list_chopped_scoreMat_list_1000bp[[1]][aamatch,]

aa_coop_dif_len_pos <- list()
aa_coop_dif_len_neg <- list()
i <- 1
# 25, 30, 35, 40, 45, 50, 60, 70, 80,
for(aa_coop_dis in c(25, 30, 35, 40, 45, 50, 60, 70, 80, 5, 10, 15, 20, 90, 100)){
  print(i)
  aa_coop_dif_len_pos[[i]] <- EnhancerChopper(Enhancer_seq = Positive_set_seq_list_char_1000bp[[1]],
                                              motifScore_output = Positive_set_motifscore_1000bp[[1]],
                                              EnhancerGR = Positive_set_gr_list_1000bp[[1]],
                                              motifList = TF.motifs.Expanded_pseudo_t,
                                              no_thresh = T,
                                              LLR_thresh = 0,
                                              piece_length = 1000,
                                              step_size = 200,
                                              normalize_by_length = F,
                                              pairwise_closeby_site_freq = T,
                                              close_thr = aa_coop_dis,
                                              LLR2PVal_lists = TF.motifs.Expanded_LLR_to_pval_pseudo,
                                              pval_threshold = 0.0001,
                                              .diff_max_LLR = F)
  aa_coop_dif_len_neg[[i]] <- EnhancerChopper(Enhancer_seq = Negative_set_seq_list_char_1000bp[[1]],
                                              motifScore_output = Negative_set_motifscore_1000bp[[1]],
                                              EnhancerGR = Negative_set_gr_list_1000bp[[1]],
                                              motifList = TF.motifs.Expanded_pseudo_t,
                                              no_thresh = T,
                                              LLR_thresh = 0,
                                              piece_length = 1000,
                                              step_size = 200,
                                              normalize_by_length = F,
                                              pairwise_closeby_site_freq = T,
                                              close_thr = aa_coop_dis,
                                              LLR2PVal_lists = TF.motifs.Expanded_LLR_to_pval_pseudo,
                                              pval_threshold = 0.0001,
                                              .diff_max_LLR = F)
  names(aa_coop_dif_len_pos)[i] <- aa_coop_dis
  names(aa_coop_dif_len_neg)[i] <- aa_coop_dis
  i <- i + 1
  
}

aa_coop_dif_len_pos <- aa_coop_dif_len_pos[sort(as.numeric(names(aa_coop_dif_len_pos)), index.return=T)$ix]
aa_coop_dif_len_neg <- aa_coop_dif_len_neg[sort(as.numeric(names(aa_coop_dif_len_neg)), index.return=T)$ix]

aa_coop_dif_len_pos_adj <- lapply(aa_coop_dif_len_pos, "[[", 4)
aa_coop_dif_len_neg_adj <- lapply(aa_coop_dif_len_neg, "[[", 4)

names(aa_coop_dif_len_pos_adj) <-names(aa_coop_dif_len_pos)
names(aa_coop_dif_len_neg_adj) <- names(aa_coop_dif_len_pos)


aa <- unlist(lapply(strsplit(names(GEMSTAT_Ensemble_train_SeqList[[2]]), split = "_"), "[[", 1))
aa_pos_t <- names(GEMSTAT_Ensemble_train_SeqList[[2]])[aa %in% "pos"]
aa_neg_t <- names(GEMSTAT_Ensemble_train_SeqList[[2]])[aa %in% "neg"]

aat_test2 <- list()
for(i in 1:length(aa_coop_dif_len_neg_adj)){
  aat_test2[[i]] <- list()
  aa_pos1 <- do.call(rbind, aa_coop_dif_len_pos_adj[[i]][aa_pos_t])
  aa_neg1 <- do.call(rbind, aa_coop_dif_len_neg_adj[[i]][aa_neg_t])
  for(j in 1:ncol(aa_pos1)){
    aa_all <- scale(c(aa_pos1[, j], aa_neg1[, j]))
    aat_test2[[i]][[j]] <- t.test(aa_all[1:nrow(aa_pos1)],
                             aa_all[(nrow(aa_pos1) + 1):length(aa_all)])
  }
}

par(mfrow = c(4, 4), mar = c(2,2,2,2))

for(i in 1:length(aat_test2)){
  plot(sort(-log10(unlist(lapply(aat_test2[[i]], "[[", 3))), decreasing = T), 
       main = names(aa_coop_dif_len_neg_adj)[i])
}



# plot -log pvalue as a function of distance for each pair of TFs
# create a matrix where each column represents a TF pair and each row represents a 
#  distance, each value is the -log10 pvlaue of the difference between nu_pairs in
#  neg vs pos training set

aa_dis_pvalmat <- matrix(nrow = length(aa_coop_dif_len_neg_adj),
                         ncol = ncol(aa_coop_dif_len_pos_adj[[1]][[1]]))
colnames(aa_dis_pvalmat) <- colnames(aa_coop_dif_len_pos_adj[[1]]$pos_1)
rownames(aa_dis_pvalmat) <- names(aa_coop_dif_len_neg_adj)
for(i in 1:nrow(aa_dis_pvalmat)){
  for(j in 1:ncol(aa_dis_pvalmat)){
    aa_dis_pvalmat[i, j] <-  -log10(aat_test2[[i]][[j]]$p.value)
  }
}
par(mfrow = c(1,1))
plot(aa_dis_pvalmat[,which(colnames(aa_dis_pvalmat) %in% "YBX1_vs_RUNX1")])

# keep column who have at least one good difference
aa_dis_pvalmat2 <- aa_dis_pvalmat > 3
aaimp <- which(colSums(aa_dis_pvalmat2) > 0)
aa_dis_pvalmat2 <- aa_dis_pvalmat[,aaimp]
aasd <- apply(aa_dis_pvalmat2, MARGIN = 2, FUN = sd)
aamean <- apply(aa_dis_pvalmat2, MARGIN = 2, FUN = mean)
aasdind <- sort(aasd, decreasing = T, index.return = T)$ix
aameanind <- sort(aamean, decreasing = T, index.return = T)$ix
aasdmean <- aasd/aamean
aasdmeanind <- sort(aasdmean, decreasing = T, index.return = T)$ix
par(mfrow = c(6, 4), mar = c(2,2,2,2))
aanames <- c("ESR1_1_vs_AR", "ESR1_1_vs_CEBPB", "FOXA1_vs_ESR1_2", "FOXA1_vs_ESR1_1", "GATA3_vs_ESR1_1", "GATA3_vs_ESR1_2",
             "JUN_1_vs_JUN_1", "JUN_2_vs_ESR1_1", "LEF1_vs_ESR1_2", "LEF1_vs_JUN_1", "MYC_vs_ESR1_1", "NKX3-1_vs_JUN_1",
             "NR5A2_vs_ESR1_1", "NR5A2_vs_ESR1_2", "PBX1_vs_ESR1_1", "PBX1_vs_ESR1_2", "PGR_vs_ESR1_1", "POU5F1_vs_ESR1_1",
             "POU5F1_vs_ESR1_2", "PPARD_vs_ESR1_1", "RUNX1_vs_ESR1_1", "YBX1_vs_ESR1_1", "YBX1_vs_JUN_1")
for(i in 1:23){
  plot(aa_dis_pvalmat2[, aanames[i]], xaxt = "n",
       main = aanames[i], type = "l", ylim =c(1,7))
  axis(side = 1, 
       at = c(1:nrow(aa_dis_pvalmat2)),
       labels = rownames(aa_dis_pvalmat2), las=2)
}



hist(unlist(lapply(aat_test2, "[[", 3)))




