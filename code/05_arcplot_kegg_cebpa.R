# Date created: 1/4/2022
# Date modified: 11/11/2022
# Author: Maren Buettner
# Project: Liver Ageing, Arc plots to visualise Gene Expression patterns for 
# certain signalling pathways in old and young CEBPA KO mice
# gene lists were curated by Kelvin Yin

# install devtools
#install.packages("devtools")
# load devtools
library(devtools)
# install arcdiagram
devtools::install_github('gastonstat/arcdiagram')
#load json file for KEGG pathways
library(rjson)

# load arcdiagram
library(arcdiagram)

library(igraph)

#genes of interest
hnf4a <- "ENSMUSG00000017950"
cebpa <- "ENSMUSG00000034957"
ctcf <-  "ENSMUSG00000005698"

#custom gene list
ref_genes <- read.table('~/Documents/Collaborations/Martinez_HPC/references/List_for Arc plot.csv',
                        sep=';' , header =TRUE
                        )


#we have computed the mean expression and co-expression of the 
#respective genes in python, then save the co-expression result as json file and
#mean gene expression as csv table.

#in addition, we read in the gene names and mean gene expression (overall)
mean_expr_all <- read.table(file = "~/Documents/Collaborations/Martinez_HPC/tables/220406_gene_mean_all.csv", 
                            header = TRUE, sep=','
)

##############
# young wt
# read result files from python.
co_corr_young <- fromJSON(file = '~/Documents/Collaborations/Martinez_HPC/tables/221111_co_expression_CEBPA KO_young.json')

mean_expr_young <- read.table(file = "~/Documents/Collaborations/Martinez_HPC/tables/221111_mean_expression_hep_CEBPA KO_young.csv", 
                              header = TRUE, sep=','
)

#filter for low mean expression
#mean_expr_young <- mean_expr_young[mean_expr_young$gene_mean_hep>1 | 
#                                     mean_expr_young$X %in% c(hnf4a, cebpa),]
mean_expr_young <- mean_expr_young[mean_expr_young$gene_mean_hep>1 | 
                                     mean_expr_young$X %in% c(hnf4a, cebpa, ctcf),]

#set min co-expression threshold
min_corr <- 0.01

#format co-expression - create an edgelist
edge_list_coex <-matrix(nrow=0, ncol=2)
edge_weights_coex <- matrix(nrow=0, ncol=1)
for (idx in 1:length(co_corr_young)){
  names_coex <- names(co_corr_young[idx])
  tmp_id <- strsplit(names_coex, '-')[[1]]
  #check whether the gene names are in the mean expression list
  if (sum(tmp_id %in% mean_expr_young$X)==2){
    num_corr <- as.numeric(co_corr_young[idx][names_coex])
    if(abs(num_corr)>=min_corr ){
      edge_list_coex <- rbind(edge_list_coex, tmp_id)
      edge_weights_coex <- rbind(edge_weights_coex, num_corr)
    }
  
  }
  
}


#filter further by the unique gene names from edgelist
edgelist_names <- unique(c(edge_list_coex[,1], edge_list_coex[,2]))
#pathway_groups <- pathway_groups[pathway_groups$genes %in% edgelist_names,]
mean_expr_young <- mean_expr_young[mean_expr_young$X %in% edgelist_names,]


edge_list_coex <- data.frame(edge_list_coex)
colnames(edge_list_coex) <- c('V1', 'V2')
edge_list_coex$V1 <- with(edge_list_coex, factor(V1, 
                                          levels = mean_expr_all$X, 
                                          labels = mean_expr_all$gene_name))
edge_list_coex$V2 <- with(edge_list_coex, factor(V2, 
                                                 levels = mean_expr_all$X, 
                                                 labels = mean_expr_all$gene_name))
#get pathway groups and filter for the unique genes that are also present in the edge_list
pathway_groups <- data.frame(ref_genes)
names(pathway_groups) <- c("genes", "group")
pathway_groups <- pathway_groups[pathway_groups$genes %in% unique(c(edge_list_coex$V1, edge_list_coex$V2)),]
#remove duplicates
pathway_groups <- pathway_groups[!duplicated(pathway_groups$genes),]


# get edgelist
edgelist <- as.matrix(edge_list_coex)
# get vertex labels (gene names)
vlabels = as.character(pathway_groups$genes)
# get vertex groups
vgroups = pathway_groups$group
# get vertex fill color
color_selection <- c("#BEBADA", "#8b91d4" ,"#67799f" ,"#d3e4a6" ,"#b2c771",
                      "#9ca243", "#657845", "#f3d581", "#f0c753","#b69f60", "#988053")
color_border <- c("#aba7c4", "#6f74a9", "#51607f", "#bdcd95" ,"#8e9f5a", "#7c8135",
                  "#5a6c3e" ,"#dabf74", "#b89e54", "#a38f56", "#745926")
vfill = as.character(factor(x=pathway_groups$group, 
               levels = unique(pathway_groups$group), 
               labels = color_selection[1:length(unique(pathway_groups$group))]
               ))
# get vertex border color
vborders = as.character(factor(x=pathway_groups$group, 
                  levels = unique(pathway_groups$group), 
                  labels = color_border[1:length(unique(pathway_groups$group))]
))
# get vertex degree
degrees = mean_expr_young$gene_mean_hep
# get edges value
values = unlist(edge_weights_coex)#get.edge.attribute(g_coex, "value")

##Order nodes
# load reshape
library(reshape)
# data frame with vgroups, degree, vlabels and ind
x = data.frame(vgroups, degrees, vlabels, vfill, vborders)
# arranging by vgroups and degrees
x = x[with(x, order(vgroups, degrees)),]
new_ord <-x$ind

# plot arc diagram
arcplot(edgelist, ordering=x$vlabels, 
        labels=vlabels, 
        cex.labels=0.8,
        show.nodes=TRUE, col.nodes=vborders, bg.nodes=vfill,
        cex.nodes = log(degrees+1), pch.nodes=21,
        lwd.nodes = 2, line=-0.5,
        col.arcs = hsv(0, 0, 0.2, 0.25), lwd.arcs = 100 * (values-min_corr))

######################3
# old wt

# read result files from python.
co_corr_old <- fromJSON(file = '~/Documents/Collaborations/Martinez_HPC/tables/221111_co_expression_CEBPA KO_old.json')

mean_expr_old <- read.table(file = "~/Documents/Collaborations/Martinez_HPC/tables/221111_mean_expression_hep_CEBPA KO_old.csv", 
                              header = TRUE, sep=','
)
#filter for the same genes as in the young 
mean_expr_old <- mean_expr_old[mean_expr_old$X %in% mean_expr_young$X,]

#set min co-expression threshold
min_corr <- 0.01

#format co-expression - create an edgelist
edge_list_coex <-matrix(nrow=0, ncol=2)
edge_weights_coex <- matrix(nrow=0, ncol=1)
for (idx in 1:length(co_corr_old)){
  names_coex <- names(co_corr_old[idx])
  tmp_id <- strsplit(names_coex, '-')[[1]]
  #check whether the gene names are in the mean expression list
  if (sum(tmp_id %in% mean_expr_old$X)==2){
    num_corr <- as.numeric(co_corr_old[idx][names_coex])
    if(abs(num_corr)>=min_corr ){
      edge_list_coex <- rbind(edge_list_coex, tmp_id)
      edge_weights_coex <- rbind(edge_weights_coex, num_corr)
    }

  }
  
}

#filter further by the unique gene names from edgelist
edgelist_names <- unique(c(edge_list_coex[,1], edge_list_coex[,2]))
#pathway_groups <- pathway_groups[pathway_groups$genes %in% edgelist_names,]
expressed_subset <- mean_expr_old$X %in% edgelist_names 
mean_expr_old <- mean_expr_old[expressed_subset,]


#reorder pathway_groups
#ord <- order(pathway_groups$genes)
#pathway_groups <- pathway_groups[ord,]

#rename gene names from ensemblID to gene symbol 
#pathway_groups$genes <- with(pathway_groups, factor(genes, 
#                           levels = mean_expr_all$X, labels = mean_expr_all$gene_name))
edge_list_coex <- data.frame(edge_list_coex)
colnames(edge_list_coex) <- c('V1', 'V2')
edge_list_coex$V1 <- with(edge_list_coex, factor(V1, 
                                                 levels = mean_expr_all$X, 
                                                 labels = mean_expr_all$gene_name))
edge_list_coex$V2 <- with(edge_list_coex, factor(V2, 
                                                 levels = mean_expr_all$X, 
                                                 labels = mean_expr_all$gene_name))
#get pathway groups and filter for the unique genes that are also present in the edge_list
pathway_groups <- data.frame(ref_genes)
names(pathway_groups) <- c("genes", "group")
pathway_groups <- pathway_groups[pathway_groups$genes %in% unique(c(edge_list_coex$V1, edge_list_coex$V2)),]
#remove duplicates
pathway_groups <- pathway_groups[!duplicated(pathway_groups$genes),]


# get edgelist
edgelist <- as.matrix(edge_list_coex)
# get vertex labels (gene names)
vlabels = as.character(pathway_groups$genes)
# get vertex groups
vgroups = pathway_groups$group
# get vertex fill color
color_selection <- c("#BEBADA", "#8b91d4" ,"#67799f" ,"#d3e4a6" ,"#b2c771",
                     "#9ca243", "#657845", "#f3d581", "#f0c753","#b69f60", "#988053")
color_border <- c("#aba7c4", "#6f74a9", "#51607f", "#bdcd95" ,"#8e9f5a", "#7c8135",
                  "#5a6c3e" ,"#dabf74", "#b89e54", "#a38f56", "#745926")
vfill = as.character(factor(x=pathway_groups$group, 
                            levels = unique(pathway_groups$group), 
                            labels = color_selection[1:length(unique(pathway_groups$group))]
))
# get vertex border color
vborders = as.character(factor(x=pathway_groups$group, 
                               levels = unique(pathway_groups$group), 
                               labels = color_border[1:length(unique(pathway_groups$group))]
))
# get vertex degree
degrees = mean_expr_old$gene_mean_hep
degrees_ref = mean_expr_young$gene_mean_hep[expressed_subset]
# get edges value
values = unlist(edge_weights_coex)#get.edge.attribute(g_coex, "value")

##Order nodes
# load reshape
library(reshape)
# data frame with vgroups, degree, vlabels and ind
x = data.frame(vgroups, degrees, vlabels, vfill, vborders, degrees_ref)
# arranging by vgroups and degrees
x = x[with(x, order(vgroups, degrees_ref)),]
new_ord <-x$ind

# plot arc diagram
arcplot(edgelist, ordering=x$vlabels, 
        labels=vlabels, 
        cex.labels=0.8,
        show.nodes=TRUE, col.nodes=vborders, bg.nodes=vfill,
        cex.nodes = log(degrees+1), pch.nodes=21,
        lwd.nodes = 2, line=-0.5,
        col.arcs = hsv(0, 0, 0.2, 0.25), lwd.arcs = 100 * (values-min_corr))

#########################
### TEST
# read 'gml' file
mis_graph = read.graph('~/Documents/R/les_miserables.gml', format="gml")
# get edgelist
edgelist = get.edgelist(mis_graph)
# get vertex labels
vlabels = get.vertex.attribute(mis_graph, "label") #Gene names
# get vertex groups
vgroups = get.vertex.attribute(mis_graph, "group") #Pathway IDs
# get vertex fill color
vfill = get.vertex.attribute(mis_graph, "fill") #colors of pathways
# get vertex border color
vborders = get.vertex.attribute(mis_graph, "border") #colors of pathways, maybe darker
# get vertex degree
degrees = degree(mis_graph) #mean expression
# get edges value
values = get.edge.attribute(mis_graph, "value") #co-expression value

##Order nodes
# load reshape
library(reshape)
# data frame with vgroups, degree, vlabels and ind
x = data.frame(vgroups, degrees, vlabels, ind=1:vcount(mis_graph))
# arranging by vgroups and degrees
y = x[order(vgroups, degrees),]
# get ordering 'ind'
new_ord = y$ind

# plot arc diagram
arcplot(edgelist, ordering=new_ord, labels=vlabels, cex.labels=0.8,
        show.nodes=TRUE, col.nodes=vborders, bg.nodes=vfill,
        cex.nodes = log(degrees)+0.5, pch.nodes=21,
        lwd.nodes = 2, line=-0.5,
        col.arcs = hsv(0, 0, 0.2, 0.25), lwd.arcs = 1.5 * values)
