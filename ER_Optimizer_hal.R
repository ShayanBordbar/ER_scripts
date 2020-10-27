#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(sigmoid)
load(Inputs.RData)
OptimizedExp <- list()
OptimizedExp[[as.integer(args[1])]] <- optim(par = runif(n = ((ncol(aaaaaa_optim_dataset$gene_conditon_mat) - 1)/2 + 1), min = -10, max = 10), fn = obj_func_linear_minEnh, my_data = aaaaaa_optim_dataset)
save.image(paste0("Optim_res_", args[1], ".RData"))