# Gene enhancer TF network construction and visualization
library(dplyr)
library(igraph)
#library(netbiov)
#devtools::install_github("frankkramer-lab/mully")
library(mully)

Enhancer_Enhancer_interac_mat
Enhancer_Enhancer_overlap_mat
Enhancer_gene_20kb_vicinity_mat
gene_gene_vicinity_mat_100kb
promoter_promoter_interaction_mat

# TF vertices
ER_KD_TF_enh_network$TF[ER_KD_TF_enh_network$TF == "ESR1_2"] <- "ESR1"
GRN_TF_vertices <- unique(ER_KD_TF_enh_network$TF)
# TF edges
TF_TF_edges <- rbind(c("ESR1", "GATA3"),
                     c("ESR1", "FOXA1"),
                     c("ESR1", "RARA"),
                     c("ESR1", "TFAP2C"),
                     c("ESR1", "PBX1"),
                     c("ESR1", "PGR"),
                     c("ESR1", "YBX1"),
                     c("RUNX1"), c("YBX1"))


# enhancer vertices
GRN_enhancer_vertices <- sort(unique(ER_KD_TF_enh_network$Enhancer))
# TF_enhacer edges: from knockdown


# TF_enhacer edges: from chip data
TF_enhancer_chip_network <- TF_enhancer_chip_network[TF_enhancer_chip_network$Enhancer %in% GRN_enhancer_vertices, ]
TF_enhancer_chip_network <- TF_enhancer_chip_network[TF_enhancer_chip_network$TF %in% GRN_TF_vertices, ]

# KD predicted interactions supported by ChIP
ER_TF_enh_KD_Chip_intersect <- inner_join(ER_KD_TF_enh_network,
                                          TF_enhancer_chip_network)

# enhancer-enhancer edges
Enhancer_Enhancer_interac_mat
Enhancer_Enhancer_overlap_mat
enhncer_enhancer_network_edges <- rbind(Enhancer_Enhancer_overlap_mat,
                                        Enhancer_Enhancer_interac_mat)
enhncer_enhancer_network_edges <- enhncer_enhancer_network_edges[!duplicated(enhncer_enhancer_network_edges),]


Enh_gene_Experimental_network <- rbind(Enhancer_gene_20kb_vicinity_mat,
                           Enhancer_genePromoter_interac_mat_experimental_only)
Enh_gene_Experimental_network <- Enh_gene_Experimental_network[!duplicated(Enh_gene_Experimental_network),]
Enh_gene_Experimental_network <- Enh_gene_Experimental_network[Enh_gene_Experimental_network$enh_name %in%GRN_enhancer_vertices, ]
# gene vertices
GRN_gene_vertices <- sort(unique(Enh_gene_Experimental_network$gene_entID))

# enhancer gene network --> overlap or interaction
Enh_gene_Experimental_network


# gene-gene nearby (100kb) , or promoter interaction
gene_gene_network_edges <- rbind(gene_gene_vicinity_mat_100kb
                                 ,promoter_promoter_interaction_mat)

gene_gene_network_edges <- gene_gene_network_edges[gene_gene_network_edges$gene1_entID %in% GRN_gene_vertices, ]
gene_gene_network_edges <- gene_gene_network_edges[gene_gene_network_edges$gene2_entID %in% GRN_gene_vertices, ]
gene_gene_network_edges <- gene_gene_network_edges[!duplicated(gene_gene_network_edges),]


GRN_vertices_list <- list(TF = GRN_TF_vertices, Enhancer = GRN_enhancer_vertices, Gene = GRN_gene_vertices)
GRN_vertices_df <- data.frame(rbind(cbind(GRN_TF_vertices, rep("TF", length(GRN_TF_vertices))), 
                         cbind(GRN_enhancer_vertices, rep("Enhancer", length(GRN_enhancer_vertices))),
                         cbind(GRN_gene_vertices, rep("Gene", length(GRN_gene_vertices)))), stringsAsFactors = F)
names(GRN_vertices_df) <- c("vertex_name", "type")
names(TF_TF_edges) <- c("V1", "V2")
names(ER_KD_TF_enh_network) <-  c("V1", "V2")
names(enhncer_enhancer_network_edges) <-  c("V1", "V2")
names(Enh_gene_Experimental_network) <-  c("V1", "V2")
names(gene_gene_network_edges) <-  c("V1", "V2")
names(ER_TF_enh_KD_Chip_intersect) <-  c("V1", "V2")

GRN_edges_final_full <- rbind(TF_TF_edges,
                              ER_KD_TF_enh_network,
                              enhncer_enhancer_network_edges,
                              Enh_gene_Experimental_network, 
                              gene_gene_network_edges
                              )
GRN_edges_final_full_type <- c(rep("TF_TF", nrow(TF_TF_edges)),
                               rep("TF_Enhancer", nrow(ER_KD_TF_enh_network)),
                               rep("Enhancer_Enhancer", nrow(enhncer_enhancer_network_edges)),
                               rep("Enhancer_Gene", nrow(Enh_gene_Experimental_network)),
                               rep("Gene_Gene", nrow(gene_gene_network_edges))
                               )
GRN_edges_final_full_type <- GRN_edges_final_full_type[GRN_edges_final_full$V1 %in% GRN_vertices_df$vertex_name]
GRN_edges_final_full <- GRN_edges_final_full[GRN_edges_final_full$V1 %in% GRN_vertices_df$vertex_name, ]
GRN_edges_final_full_type <- GRN_edges_final_full_type[GRN_edges_final_full$V2 %in% GRN_vertices_df$vertex_name]
GRN_edges_final_full <- GRN_edges_final_full[GRN_edges_final_full$V2 %in% GRN_vertices_df$vertex_name, ]


GRN_edges_final_ChiPSUPP <- rbind(TF_TF_edges,
                                  ER_TF_enh_KD_Chip_intersect, 
                                  #enhncer_enhancer_network_edges, 
                                  Enh_gene_Experimental_network#, 
                                  #gene_gene_network_edges
                                  )
GRN_edges_final_ChiPSUPP_type <- c(rep("TF_TF", nrow(TF_TF_edges)),
                               rep("TF_Enhancer", nrow(ER_TF_enh_KD_Chip_intersect)),
                              # rep("Enhancer_Enhancer", nrow(enhncer_enhancer_network_edges)),
                               rep("Enhancer_Gene", nrow(Enh_gene_Experimental_network))#,
                               #rep("Gene_Gene", nrow(gene_gene_network_edges))
                              )

GRN_edges_final_ChiPSUPP <- GRN_edges_final_ChiPSUPP[GRN_edges_final_ChiPSUPP$V1 %in% GRN_vertices_df$vertex_name, ]
GRN_edges_final_ChiPSUPP_type <- GRN_edges_final_ChiPSUPP_type[GRN_edges_final_ChiPSUPP$V1 %in% GRN_vertices_df$vertex_name]
GRN_edges_final_ChiPSUPP <- GRN_edges_final_ChiPSUPP[GRN_edges_final_ChiPSUPP$V2 %in% GRN_vertices_df$vertex_name, ]
GRN_edges_final_ChiPSUPP_type <- GRN_edges_final_ChiPSUPP_type[GRN_edges_final_ChiPSUPP$V2 %in% GRN_vertices_df$vertex_name]



# TF_Enh_gene_network_KD <- full_join(x = ER_KD_TF_enh_network,
#                             y = Enh_gene_Experimental_network,
#                             by = c("Enhancer"="enh_name"))
# TF_Enh_gene_network_KD_filtered_gene <- TF_Enh_gene_network_KD[!is.na(TF_Enh_gene_network_KD$gene_entID),]
# 
# 
# TF_Enh_gene_network_KD_Chip_intersect <- left_join(x = ER_TF_enh_KD_Chip_intersect,
#                                     y = Enh_gene_Experimental_network,
#                                     by = c("Enhancer"="enh_name"))
# TF_Enh_gene_network_KD_Chip_intersect_filteredGene <- TF_Enh_gene_network_KD_Chip_intersect[!is.na(TF_Enh_gene_network_KD_Chip_intersect$gene_entID),]
# 
# 
# 
# TF_Enh_gene_network_KD_filtered_gene
# TF_Enh_gene_network_KD_Chip_intersect_filteredGene



#aagrph <- graph_from_data_frame(d =GRN_edges_final_full, directed = F, vertices = GRN_vertices_df)
data("color_list")


save(list = c("GRN_edges_final_full", 
              "GRN_edges_final_ChiPSUPP", 
              "GRN_vertices_list",
              "GRN_vertices_df",
              "TF_TF_edges",
              "ER_TF_enh_KD_Chip_intersect", 
              "ER_KD_TF_enh_network",
              "enhncer_enhancer_network_edges", 
              "gene_gene_network_edges", 
              "GRN_edges_final_ChiPSUPP_type", 
              "GRN_edges_final_full_type"), 
     file = "~/Documents/Shayan/BioInf/EstrogenReceptor/GRN_workplace.RData")




# plotting using netbiov

aacol <- RColorBrewer::brewer.pal(n = 12, name = "Set3")
#color for edges
aaedcol <- character(length(GRN_edges_final_full_type))
aaaaun <- unique(GRN_edges_final_full_type)
for(i in 1:length(aaaaun)){
  aaedcol[GRN_edges_final_full_type == aaaaun[i]] <- aacol[i]
}

#hc <- rgb(t(col2rgb(heat.colors(20)))/255,alpha=.2) 
aacl <- character(nrow(GRN_vertices_df))
aasize <- numeric(nrow(GRN_vertices_df))

aacl[GRN_vertices_df$type == "TF"] <- aacol[5]
aacl[GRN_vertices_df$type == "Enhancer"] <- aacol[6]
aacl[GRN_vertices_df$type == "Gene"] <- aacol[7]

plot(c(1:10),pch = 16, cex = 2, col = c(aacol[5],  aacol[6],  aacol[7]) )
aatable <- table(GRN_vertices_df$type)

aasize[GRN_vertices_df$type == "TF"] <- 0.5
aasize[GRN_vertices_df$type == "Enhancer"] <- 0.1
aasize[GRN_vertices_df$type == "Gene"] <- 0.1




#drawing a sample of this network with 100  edges
#aa_Samp <- sample(c(1:nrow(GRN_edges_final_full)), 100)
#aagrph <- graph_from_data_frame(d =GRN_edges_final_full, directed = F, vertices = GRN_vertices_df)
#aaedcol <- aaedcol[aa_Samp]
png(filename = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Plots/GRN_Test1.png", 
    width = 10*300, 
    height = 10*300, 
    res = 300,
    pointsize = 10)
xx <- mst.plot.mod(aagrph, 
                   #vertex.color= aacl,
                   v.size=aasize,
                   sf=-20, 
                   colors=aaedcol,
                   e.size=.5,
                   mst.e.size=.75,
                   e.lab = F,
                   layout.function=layout.fruchterman.reingold)


dev.off()







# draw modular one
GRN_edges_final_full_noGGEE <- GRN_edges_final_full[! GRN_edges_final_full_type %in% c("Enhancer_Enhancer", "Gene_Gene"),]
aagrph <- graph_from_edgelist(el =  as.matrix(GRN_edges_final_full_noGGEE),
                                directed = F)
GRN_edges_final_full_type_filtered <- GRN_edges_final_full_type[! GRN_edges_final_full_type %in% c("Enhancer_Enhancer", "Gene_Gene")]
aaedcol <- character(length(GRN_edges_final_full_type_filtered))
aaaaun <- unique(GRN_edges_final_full_type_filtered)
for(i in 1:length(aaaaun)){
  aaedcol[GRN_edges_final_full_type_filtered == aaaaun[i]] <- color.list$morning1[i]
}

aavname <- V(aagrph)$name
aavname_type <- GRN_vertices_df$type[match(aavname, GRN_vertices_df$vertex_name)]
aacl <- character(length = (aavname))
aasize <- numeric(length = (aavname))

aacl[aavname_type == "TF"] <- color.list$warm[1]
aacl[aavname_type == "Enhancer"] <- color.list$bluebush[1]
aacl[aavname_type == "Gene"] <- color.list$lily[10]
aasize[aavname_type == "TF"] <- 0.5
aasize[aavname_type == "Enhancer"] <- 0.1
aasize[aavname_type == "Gene"] <- 0.1


png(filename = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Plots/GRN_Test3.png", 
    width = 10*300, 
    height = 10*300, 
    res = 300,
    pointsize = 10)
xx2 <- plot.modules(aagrph, 
                    layout.function = c(layout.fruchterman.reingold),
                    modules.color = sample(color.list$bright), 
                    #layout.overall = layout.star ,
                   vertex.color=aacl,
                   v.size=aasize,
                   sf=40,
                   tkplot=FALSE,
                   colors=aaedcol,
                   e.size=.5,
                   mst.e.size=.75,
                   e.lab = F)

dev.off()


#data("PPI_Athalina")
data("color_list")
id <- plot.modules(aagrph, layout.function = c(layout.fruchterman.reingold),
                   modules.color = sample(color.list$bright),
                   layout.overall = layout.star ,
                   sf=40, tkplot=FALSE)
# draw a hirarchical one

png(filename = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Plots/GRN_Test4.png", 
    width = 10*300, 
    height = 10*300, 
    res = 300,
    pointsize = 10)
id <- plot.modules(aagrph, 
                   layout.function = layout.kamada.kawai ,
                   sf = 0, 
                   col.grad=list(color.list$citynight),
                   tkplot=FALSE,
                   vertex.color=aacl,
                   v.size=aasize,
                   colors=aaedcol,
                   e.size=.5)
dev.off()

# draw for one TF at a time

png(filename = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Plots/GRN_Test5.png", 
    width = 10*300, 
    height = 10*300, 
    res = 300,
    pointsize = 10)
xx <- plot.spiral.graph(aagrph, tp = 179,vertex.color=aacl )
dev.off()



png(filename = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Plots/GRN_Test5.png", 
    width = 10*300, 
    height = 10*300, 
    res = 300,
    pointsize = 10)
xx <- plot.abstract.module(aagrph, tkplot = FALSE, layout.function = layout.star,vertex.color=aacl)
dev.off()

?level.plot

aavname <- V(aagrph)$name
aavname_type <- GRN_vertices_df$type[match(aavname, GRN_vertices_df$vertex_name)]
aacl <- character(length = (aavname))
aasize <- numeric(length = (aavname))

aavname_TFs <- aavname[aavname_type == "TF"]
as.numeric(V(aagrph)[aavname_TFs])
aacl[aavname_type == "TF"] <- color.list$warm[1]
aacl[aavname_type == "Enhancer"] <- color.list$bluebush[1]
aacl[aavname_type == "Gene"] <- color.list$lily[10]
aasize[aavname_type == "TF"] <- 0.5
aasize[aavname_type == "Enhancer"] <- 0.1
aasize[aavname_type == "Gene"] <- 0.1


png(filename = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Plots/GRN_Test6.png", 
    width = 10*300, 
    height = 10*300, 
    res = 300,
    pointsize = 10)
xx <- level.plot(aagrph, 
                 tkplot=FALSE, 
                 level.spread=FALSE, 
                 layout.function=layout.fruchterman.reingold,
                 type =1,
                 #initial_nodes = which(GRN_vertices_df$type == "TF"), 
                 nodeset= list(which(aavname_type == "TF"),
                               which(aavname_type == "Gene"))
                 #,vertex.color=aacl
                 )
dev.off()
layout.kamada.kawai
                           

d



png(filename = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Plots/GRN_Test7.png", 
    width = 10*300, 
    height = 10*300, 
    res = 300,
    pointsize = 10)

aavid <- (V(aagrph)[aa_wtf])
#aavid <- as.numeric(V(aagrph)$name)
aa_wtf <- which(aavname_type == "TF") 
aa_wge <- which(aavname_type == "Gene")
plot.igraph(aagrph)
id <- level.plot(x = aagrph, 
                 #layout.function = layout.fruchterman.reingold,
                 #type = 1, 
                 #col.s1 = color.list$bluebush[4],
                 #col.s2 = "pink", 
                 #v.size = 0.5,
                 #e.size = 0.5, 
                # v.lab = T, 
                 #v.lab.cex = 0.3,
                #init_nodes=16,
                 initial_nodes= c(1,2,3,4,5,6) + 2#,
                 #vertex.colors=c("white", "white", "white"),
                # edge.col=c("grey", "grey", "grey", "grey"),
                # nodeset= list(aa_wtf,aa_wge),
                 #tkplot=FALSE,
                #level.spread=TRUE,
                # order_degree="out"
                )

dev.off()
##########################################################################################
# using "mully"

# g <- mully("MyFirstMully",direct = F)
# g <- addLayer(g, c("gene", "disease", "drug"))
# 
# g=addNode(g,"d1","disease",attributes=list(type="t1"))
# print("Node d1 added as disease")
# 
# g=addNode(g,"d2","disease",attributes=list(type="t1"))
# print("Node d2 added as disease")
# 
# g=addNode(g,"d3","disease",attributes=list(type="t1"))
# print("Node d3 added as disease")
# 
# g=addNode(g,"dr1","drug",attributes=list(effect="strong"))
# print("Node dr1 added as drug")
# 
# g=addNode(g,"dr2","drug",attributes=list(effect="strong"))
# print("Node dr2 added as drug")
# 
# g=addNode(g,"dr3","drug",attributes=list(effect="moderate"))
# print("Node dr3 added as drug")
# 
# g=addNode(g,"g1","gene",attributes=list(desc="AF"))
# print("Node g1 added as gene")
# 
# g=addNode(g,"g2","gene",attributes=list(desc="BE"))
# print("Node g2 added as gene")
# 
# #See vertices attributes
# print(getNodeAttributes(g))
# 
# g=addEdge(g,"dr1","d2",list(name="treats"))
# g=addEdge(g,"dr1","d2",list(name="extraEdge"))
# g=addEdge(g,"d2","g1",list(name="targets"))
# g=addEdge(g,"g2","dr3",list(name="mutates and causes"))
# g=addEdge(g,"dr3","d3",list(name="treats"))
# 
# print(getEdgeAttributes(g))
# removeEdge(g,"d2","dr1",multi=T)
# 

#Merge both graphs
g12=merge(g,g1)

plot(g,layout = "scaled")
plot3d(g)


# GRN_vertices_df$type[GRN_vertices_df$type == "Enhancer"] <- "enhancer"
# GRN_vertices_df$type[GRN_vertices_df$type == "TF"] <- "tf"
# GRN_vertices_df$type[GRN_vertices_df$type == "Gene"] <- "gene"

# GRN_edges_final_full$V1 <- levels(GRN_edges_final_full$V1 )[as.numeric(GRN_edges_final_full$V1)]
# GRN_edges_final_full$V2 <- levels(GRN_edges_final_full$V2 )[as.numeric(GRN_edges_final_full$V2)]

aacol <- RColorBrewer::brewer.pal(n = 12, name = "Set3")
#color for edges
aaedcol <- character(length(GRN_edges_final_full_type))
aaaaun <- unique(GRN_edges_final_full_type)
for(i in 1:length(aaaaun)){
  aaedcol[GRN_edges_final_full_type == aaaaun[i]] <- aacol[i+ 3]
}
plot(c(1:10),pch = 16, cex = 2, col = c(aacol[5],  aacol[6],  aacol[7]) )
#hc <- rgb(t(col2rgb(heat.colors(20)))/255,alpha=.2)

#hc <- rgb(t(col2rgb(heat.colors(20)))/255,alpha=.2) 
aacl <- character(nrow(GRN_vertices_df))
aasize <- numeric(nrow(GRN_vertices_df))

aacl[GRN_vertices_df$type == "TF"] <- aacol[5]
aacl[GRN_vertices_df$type == "Enhancer"] <- aacol[6]
aacl[GRN_vertices_df$type == "Gene"] <- aacol[7]

plot(c(1:10),pch = 16, cex = 2, col = c(aacol[5],  aacol[6],  aacol[7]) )
aatable <- table(GRN_vertices_df$type)

aasize[GRN_vertices_df$type == "TF"] <- 10
aasize[GRN_vertices_df$type == "Enhancer"] <- 0.1
aasize[GRN_vertices_df$type == "Gene"] <- 0.1


GRN_graph <- mully(name = "ER_GRN", direct = T)
GRN_graph <- addLayer(GRN_graph, c("tf", "enhancer", "gene"))

for(aa_cur_node in 1:nrow(GRN_vertices_df)){
  print(aa_cur_node)
  GRN_graph=addNode(GRN_graph,
                    GRN_vertices_df$vertex_name[aa_cur_node],
                    GRN_vertices_df$type[aa_cur_node],
                    attributes=list(size=aasize[aa_cur_node],
                                    color = aacl[aa_cur_node]))
}
for(aa_cir_edge in 1:nrow(GRN_edges_final_full)){
  print(aa_cir_edge)
  GRN_graph=addEdge(GRN_graph,
                    GRN_edges_final_full$V1[aa_cir_edge],
                    GRN_edges_final_full$V2[aa_cir_edge],
                    list(name=GRN_edges_final_full_type[aa_cir_edge],
                         color = aaedcol[aa_cir_edge]))
}



png(filename = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Plots/GRN_Test10.png", 
    width = 10*300, 
    height = 10*300, 
    res = 300,
    pointsize = 10)
plot(GRN_graph,layout = "scaled")


dev.off()

#############
#Visualizing a tfap2c coop network

TFAP2C_coop_network$TF[TFAP2C_coop_network$TF == "ESR1_2"] <- "ESR1"

TFAP2C_coop_network_edges <- GRN_edges_final_full[(GRN_edges_final_full$V2 %in% TFAP2C_coop_network$Enhancer)
                                                  | (GRN_edges_final_full$V1 %in% TFAP2C_coop_network$Enhancer),]
TFAP2C_coop_network_edges_type <- GRN_edges_final_full_type[(GRN_edges_final_full$V2 %in% TFAP2C_coop_network$Enhancer)
                                                       | (GRN_edges_final_full$V1 %in% TFAP2C_coop_network$Enhancer)]
TFAP2C_coop_network_edges <- rbind(GRN_edges_final_full[c(2,4, 7),], TFAP2C_coop_network_edges)
TFAP2C_coop_network_edges_type <- c(GRN_edges_final_full_type[c(2,4, 7)], TFAP2C_coop_network_edges_type)
TFAP2C_coop_network_edges$V1 <- levels(TFAP2C_coop_network_edges$V1)[as.numeric(TFAP2C_coop_network_edges$V1)]
TFAP2C_coop_network_edges$V2 <- levels(TFAP2C_coop_network_edges$V2)[as.numeric(TFAP2C_coop_network_edges$V2)]

TFAP2C_GRN_vertices_df <- GRN_vertices_df[((GRN_vertices_df$vertex_name %in% TFAP2C_coop_network_edges$V1) |
                                             (GRN_vertices_df$vertex_name %in% TFAP2C_coop_network_edges$V2)),]
#GRN_edges_final_full[c(2,4, 7),]

TFAP2C_GRN_vertices_df$type[TFAP2C_GRN_vertices_df$type == "Enhancer"] <- "enhancer"
TFAP2C_GRN_vertices_df$type[TFAP2C_GRN_vertices_df$type == "TF"] <- "tf"
TFAP2C_GRN_vertices_df$type[TFAP2C_GRN_vertices_df$type == "Gene"] <- "gene"

# GRN_edges_final_full$V1 <- levels(GRN_edges_final_full$V1 )[as.numeric(GRN_edges_final_full$V1)]
# GRN_edges_final_full$V2 <- levels(GRN_edges_final_full$V2 )[as.numeric(GRN_edges_final_full$V2)]

#color for edges
aaedcol <- character(length(TFAP2C_coop_network_edges_type))
aaaaun <- unique(TFAP2C_coop_network_edges_type)
for(i in 1:length(aaaaun)){
  aaedcol[TFAP2C_coop_network_edges_type == aaaaun[i]] <- aacol[i+3]
}

#hc <- rgb(t(col2rgb(heat.colors(20)))/255,alpha=.2) 
aacl <- character(nrow(TFAP2C_GRN_vertices_df))
aasize <- numeric(nrow(TFAP2C_GRN_vertices_df))

aacl[TFAP2C_GRN_vertices_df$type == "tf"] <- aacol[5]
aacl[TFAP2C_GRN_vertices_df$type == "enhancer"] <- aacol[6]
aacl[TFAP2C_GRN_vertices_df$type == "gene"] <- aacol[7]

plot(c(1:10),pch = 16, cex = 2, col = c(color.list$bright[110],  color.list$bright[140],  color.list$bright[20]) )
aatable <- table(TFAP2C_GRN_vertices_df$type)

aasize[TFAP2C_GRN_vertices_df$type == "tf"] <- 10
aasize[TFAP2C_GRN_vertices_df$type == "enhancer"] <- 5
aasize[TFAP2C_GRN_vertices_df$type == "gene"] <- 5

aagnsymbol <- select(x = org.Hs.eg.db,
               keys = TFAP2C_GRN_vertices_df$vertex_name[TFAP2C_GRN_vertices_df$type == "gene"],
               columns = c("SYMBOL"),
               keytype = "ENTREZID")

all(aagnsymbol$ENTREZID == TFAP2C_GRN_vertices_df$vertex_name[TFAP2C_GRN_vertices_df$type == "gene"])


TFAP2C_GRN_vertices_df$vertex_name[TFAP2C_GRN_vertices_df$type == "gene"] <- aagnsymbol$SYMBOL

for(i in 1:nrow(aagnsymbol)){
  TFAP2C_coop_network_edges$V2[TFAP2C_coop_network_edges$V2 == aagnsymbol$ENTREZID[i]] <- aagnsymbol$SYMBOL[i]
}

TFAP2C_coop_graph <- mully(name = "TFAP2C_GRN", direct = T)
TFAP2C_coop_graph <- addLayer(TFAP2C_coop_graph, c("tf", "enhancer", "gene"))

for(aa_cur_node in 1:nrow(TFAP2C_GRN_vertices_df)){
  TFAP2C_coop_graph=addNode(TFAP2C_coop_graph,
                    TFAP2C_GRN_vertices_df$vertex_name[aa_cur_node],
                    TFAP2C_GRN_vertices_df$type[aa_cur_node],
                    attributes=list(size=aasize[aa_cur_node],
                                    color = aacl[aa_cur_node]))
}
for(aa_cir_edge in 1:nrow(TFAP2C_coop_network_edges)){
  TFAP2C_coop_graph=addEdge(TFAP2C_coop_graph,
                    TFAP2C_coop_network_edges$V1[aa_cir_edge],
                    TFAP2C_coop_network_edges$V2[aa_cir_edge],
                    list(name=TFAP2C_coop_network_edges_type[aa_cir_edge],
                         color = aaedcol[aa_cir_edge]))
}



png(filename = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Plots/GRN_Test_TFAP2C2.png", 
    width = 10*300, 
    height = 10*300, 
    res = 300,
    pointsize = 10)
plot(TFAP2C_coop_graph,layout = "scaled")


dev.off()

plot3d(TFAP2C_coop_graph, layers = F, vertex.label = TFAP2C_GRN_vertices_df$vertex_name, edge.arrow.size = 2)#, 
     #  vertex.label.color = )




save(list = c("TFAP2C_coop_graph", 
              "TFAP2C_GRN_vertices_df", 
              "TFAP2C_coop_network_edges"), 
     file = "~/Documents/Shayan/BioInf/EstrogenReceptor/GRN_TFAP2C_workplace.RData")

