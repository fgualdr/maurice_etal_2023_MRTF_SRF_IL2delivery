library( openxlsx )
library(eulerr)
library(viridis)
library(AnnotationDbi)
library(fitdistrplus)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(FactoMineR)
library(uwot)
library(FNN)
library(igraph)
library(purrr)

pval = 0.05
l2fc = log2(2)

HG_fun_upper <- function(selected,set,universe){
    # given thre vectors of enhancers names (e.g. coordinates in the form chr:start-end)
    # will compute the HG test
    if(length(selected)!=0 & length(set)!=0){
        q=length(intersect(selected,set)) 
        m=length(set) 
        n=length(universe)-m 
        k=length(selected) 
        return(phyper(         q = q, ## number of thigs you identified which are "successes" e.g. intersection
                                m = m, ## maximum possible number of successes e.g. n° of ATLAS CHIP positive
                                n = n, ## stuff wot ain't a success e.g. all the enhancers not bound by the ATLAS CHIP
                                k = k, ## enhancers affected by the inhibitor (divided into UP or DOWN)
                                lower.tail = FALSE) ) 
    }else{
        return(1)
    }
}

HG_fun_lower <- function(selected,set,universe){
    # given thre vectors of enhancers names (e.g. coordinates in the form chr:start-end)
    # will compute the HG test
    if(length(selected)!=0 & length(set)!=0){
        q=length(intersect(selected,set)) 
        m=length(set) 
        n=length(universe)-m 
        k=length(selected) 
        return(phyper(         q = q, ## number of thigs you identified which are "successes" e.g. intersection
                                m = m, ## maximum possible number of successes e.g. n° of ATLAS CHIP positive
                                n = n, ## stuff wot ain't a success e.g. all the enhancers not bound by the ATLAS CHIP
                                k = k, ## enhancers affected by the inhibitor (divided into UP or DOWN)
                                lower.tail = TRUE) ) 
    }else{
        return(1)
    }
}

JI_fun <- function(l1,l2){
    # Jaccard computation
    inter = length(intersect(l1,l2))
    #un = length((ki_enh)) #
    un = length(union(l1,l2))
    return(inter/un)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

rows_in_lim <- function(m,ld,lu){
        m[m<ld] = 0
        m[m>lu] = 0
        m[m!=0] = 1
        return(rownames(m)[which(rowSums(m) == ncol(m))])
}
split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])


sigmoid2x_impulse <- function(x, PP_param,t) {
    return(
    	(1/PP_param[3]) *
            (PP_param[2] + (PP_param[3] - PP_param[2]) *
                 (1/(1 + exp(-PP_param[1] * (t - PP_param[5]))))) *
            (PP_param[4] + (PP_param[3] - PP_param[4]) *
                 (1/(1 + exp(PP_param[1] * (t - PP_param[6])))))
    	)
}

sigmoid_monot <- function(x, PP_param,t) {
	return(
        PP_param[2] + (PP_param[3] - PP_param[2]) *
            (1/(1 + exp(-PP_param[1] * (t - PP_param[4]))))
            )
    }

Constant <- function(x, PP_param) {
	return(
        rep(PP_param[1],length(x))
            )
    }

indiv_list = function(l1,l2){
    l1l2 = intersect(l1,l2)
    l1u = l1[ l1 %!in% l1l2]
    l2u = l2[ l2 %!in% l1l2]
    return(list(l1u=l1u,l2u=l2u,l1l2=l1l2))
}

###############

sqlite_path = "MM10_Txdb_mergedTx.sqlite"
Ref_Txdb <- loadDb( sqlite_path )

anno_path = "MM10_merged_annotation.txt"
Ref_Anno <- read.delim( anno_path,sep="\t" ,row.names=1)
Ref_Anno$id = paste0(rownames(Ref_Anno),"_",Ref_Anno$Entrez.Gene.ID)
Ref_Anno_rn = Ref_Anno
rownames(Ref_Anno_rn) = Ref_Anno_rn$id

###############
msig = "/msigdb_v7.5.1_b/MSIG_list_symbol_entrez_mouse"
msig <- load( msig ) # MSIG_list_symbol_entrez_mouse

MSIG_list_symbol_entrez_mouse = lapply(MSIG_list_symbol_entrez_mouse,function(x){
    return(lapply(x,function(y){
        return(y[which(y %in% rownames(Ref_Anno_rn))])
    }))
    
})

Universe = unique(unlist(MSIG_list_symbol_entrez_mouse))

nm = names(MSIG_list_symbol_entrez_mouse)
nm = nm[grep("KEGG|ImmuneSigDB|CC: subset of GO|MF: subset of GO|BIOCARTA|BP: subset of GO|hallmark|positional|curated|chemical and genetic perturbations|PID subset of CP|Canonical pathways|WikiPathways|^C3|^MIR|^miRDB|^MIR_Legacy|^TFT|^TFT_Legacy|^GTRD|^C4|^CGN|^CM|^HPO|^C6|^C7|^VAX|^C8",nm)]

MSIG_list_symbol_entrez_mouse = MSIG_list_symbol_entrez_mouse[which( names(MSIG_list_symbol_entrez_mouse) %!in% nm)]


###############
# Import Cyril data
# Import Gualdrini et al. data
Cyril = read.xlsx("/Users/ieo5244/Documents/RT_DATA_MRTF/Supplemental_Tables_1-8Esnault2014.xlsx", "Table S2",rowNames=TRUE,colNames=TRUE,startRow=9)
Gualdrini = read.xlsx("/Users/ieo5244/Documents/RT_DATA_MRTF/Supplemental_Tables_1Gualdrini2016.xlsx", "Table S1 RNAseq",rowNames=TRUE,colNames=TRUE,startRow=2)

Cyril_list = list(
    FCS_induced = Cyril$Name_conversion[which(Cyril[,"Induced.by.FCS"] == 1 )],
    FCS_induced_sensitive_linked_signals = Cyril$Name_conversion[which(Cyril[,"Induced.by.FCS.and.sensitive.to.SRF-linked.signal"] == 1 )],
    Inhibited_Lat_pm_UO = Cyril$Name_conversion[which(Cyril[,"Inhibited.by.LatB.±.U0126"] == 1 )],
    Inhibited_Lat_pm_UO_direct = Cyril$Name_conversion[which(Cyril[,"Inhibited.by.LatB.±.U0126"] == 1 & Cyril[,"Direct:.within.gene.feature.or.<2kb.from.TSS"] == 1 )],
    Inhibited_UO_pm_Lat = Cyril$Name_conversion[which(Cyril[,"Inhibited.by.U0126.±.LatB"] == 1 )],
    Inhibited_UO_pm_Lat_direct = Cyril$Name_conversion[which(Cyril[,"Inhibited.by.U0126.±.LatB"] == 1 & Cyril[,"Direct:.within.gene.feature.or.<2kb.from.TSS"] == 1 )]
)

Gualdrini_list = list(
    TPA_induced = Gualdrini$Name_conversion[which(Gualdrini[,"TPA.INDUCED"] == 1 )],
    TPA_induced_TCF_dep = Gualdrini$Name_conversion[which(Gualdrini[,"TCF.DEPEDNENT"] == 1 )],
    TPA_induced_near_or_linked = Gualdrini$Name_conversion[which(Gualdrini[,"TPA-induced,.near.or.linked.to.SRF"] == 1 )],
    TPA_induced_TCF_dep_direct = Gualdrini$Name_conversion[which(Gualdrini[,"TPA-induced,.TCF-dependent,.DIRECT.-.near.or.linked.to.SRF"] == 1 )],
    BL_enhanced_near_or_linked = Gualdrini$Name_conversion[which(Gualdrini[,"Baseline.enhanced.in.TKO,.near.or.linked.to.SRF"] == 1 )],
    BL_downreg_near_or_linked = Gualdrini$Name_conversion[which(Gualdrini[,"Baseline.reduced.in.TKO,.near.or.linked.to.SRF"] == 1 )]
)

# MSIG_list_symbol_entrez_mouse[["Cyril_list"]] = Cyril_list
# MSIG_list_symbol_entrez_mouse[["Gualdrini_list"]] = Gualdrini_list

Cyril_FCSinduced_sensitiveSRFlinked =Cyril$Name_conversion[which(Cyril[,"Inhibited.by.LatB.±.U0126"] == 1 & Cyril[,"Direct:.within.gene.feature.or.<2kb.from.TSS"] == 1 )]
Gualdrini_SRFdirect = Gualdrini$Name_conversion[which(Gualdrini[,"TPA-induced,.TCF-dependent,.DIRECT.-.near.or.linked.to.SRF"] == 1 )]

### diane : 
out_dir = "./maurice_etal_2023_MRTF_SRF_IL2delivery/"
deg_design = read.delim("./maurice_etal_2023_MRTF_SRF_IL2delivery/DEG_DESIGN",sep="\t")
diane_excel = "./maurice_etal_2023_MRTF_SRF_IL2delivery/quantifcation_exon_intron_mm10_ucsc_igenomes_diane.maurice_asf-RN22062.xlsx"
Norm_mat = read.delim("./maurice_etal_2023_MRTF_SRF_IL2delivery/Normalisation_Parameters_DIANE_M.txt",sep="\t",row.names=1)
# Make individual data_tables all, intron, exons

all = read.xlsx(diane_excel, "all",rowNames=TRUE,colNames=TRUE)
exon = read.xlsx(diane_excel, "exon",rowNames=TRUE,colNames=TRUE)
intron = all - exon[rownames(all),colnames(all)]

un = unique(stringr::str_count(deg_design$Sample_Condition,"_")) + 1
condition = stringr::str_split_fixed(deg_design$Sample_Condition, "_", un)
colnames(condition) = paste0("Condition_",1:ncol(condition))

deg_design = as.data.frame(cbind(deg_design,condition))
rownames(deg_design) = deg_design$Sample_ID

list_mat = list(norm_all=all,norm_exon=exon,norm_intron=intron)

load(paste0(out_dir,"/ListDat.RData")) ## this load resALL - list of results: we consider "all" and "intron" reads:
resALL = resALL[c("norm_all","norm_intron")]

#########
# DEG selections:
out_dir_plots = paste0(out_dir,"Plots")
dir.create(out_dir_plots)

#### Naive TCR24 RestedIL12 #########
# Extract the differences ###########
#####################################
# extract changes that are transition specific:

WT_naive_to_TCR24 = list()
up_all_WT_naive_to_TCR24 <- rownames(resALL[["norm_all"]][["WT_TCR24h_vs_WT_naive"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_all_WT_naive_to_TCR24 <-  rownames(resALL[["norm_all"]][["WT_TCR24h_vs_WT_naive"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_WT_naive_to_TCR24 <- rownames(resALL[["norm_intron"]][["WT_TCR24h_vs_WT_naive"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_WT_naive_to_TCR24 <-  rownames(resALL[["norm_intron"]][["WT_TCR24h_vs_WT_naive"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
WT_naive_to_TCR24[["up"]] = union(up_all_WT_naive_to_TCR24,up_intron_WT_naive_to_TCR24)
WT_naive_to_TCR24[["down"]] = union(down_all_WT_naive_to_TCR24,down_intron_WT_naive_to_TCR24)

KO_naive_to_TCR24 = list()
up_all_KO_naive_to_TCR24 <- rownames(resALL[["norm_all"]][["KO_TCR24h_vs_KO_naive"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_all_KO_naive_to_TCR24 <-  rownames(resALL[["norm_all"]][["KO_TCR24h_vs_KO_naive"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_KO_naive_to_TCR24 <- rownames(resALL[["norm_intron"]][["KO_TCR24h_vs_KO_naive"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_KO_naive_to_TCR24 <-  rownames(resALL[["norm_intron"]][["KO_TCR24h_vs_KO_naive"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
KO_naive_to_TCR24[["up"]] = union(up_all_KO_naive_to_TCR24,up_intron_KO_naive_to_TCR24)
KO_naive_to_TCR24[["down"]] = union(down_all_KO_naive_to_TCR24,down_intron_KO_naive_to_TCR24)

s4 <- list( "up_WT_naive_to_TCR24" = WT_naive_to_TCR24[["up"]] ,
            "up_KO_naive_to_TCR24" = KO_naive_to_TCR24[["up"]]
            )

pdf(paste0(out_dir_plots,"/01_nive_tcr_transition_up_in_WT_KO.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

s4 <- list( "down_WT_naive_to_TCR24" = WT_naive_to_TCR24[["down"]] ,
            "down_KO_naive_to_TCR24" = KO_naive_to_TCR24[["down"]]
            )

pdf(paste0(out_dir_plots,"/01_nive_tcr_transition_down_in_WT_KO.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

WT_naive_vs_KO_naive = list()
up_all_WT_naive_vs_KO_naive <- rownames(resALL[["norm_all"]][["KO_naive_vs_WT_naive"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_all_WT_naive_vs_KO_naive <-  rownames(resALL[["norm_all"]][["KO_naive_vs_WT_naive"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_WT_naive_vs_KO_naive <- rownames(resALL[["norm_intron"]][["KO_naive_vs_WT_naive"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_WT_naive_vs_KO_naive <-  rownames(resALL[["norm_intron"]][["KO_naive_vs_WT_naive"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
WT_naive_vs_KO_naive[["up"]] = union(up_all_WT_naive_vs_KO_naive,up_intron_WT_naive_vs_KO_naive)
WT_naive_vs_KO_naive[["down"]] = union(down_all_WT_naive_vs_KO_naive,down_intron_WT_naive_vs_KO_naive)

WT_tcr24_vs_KO_tcr24 = list()
up_all_WT_tcr24_vs_KO_tcr24 <- rownames(resALL[["norm_all"]][["KO_TCR24h_vs_WT_TCR24h"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_all_WT_tcr24_vs_KO_tcr24 <-  rownames(resALL[["norm_all"]][["KO_TCR24h_vs_WT_TCR24h"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_WT_tcr24_vs_KO_tcr24 <- rownames(resALL[["norm_intron"]][["KO_TCR24h_vs_WT_TCR24h"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_WT_tcr24_vs_KO_tcr24 <-  rownames(resALL[["norm_intron"]][["KO_TCR24h_vs_WT_TCR24h"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
WT_tcr24_vs_KO_tcr24[["up"]] = union(up_all_WT_tcr24_vs_KO_tcr24,up_intron_WT_tcr24_vs_KO_tcr24)
WT_tcr24_vs_KO_tcr24[["down"]] = union(down_all_WT_tcr24_vs_KO_tcr24,down_intron_WT_tcr24_vs_KO_tcr24)

s4 <- list( "up_WT_naive_vs_KO_naive" = WT_naive_vs_KO_naive[["up"]] ,
            "up_WT_tcr24_vs_KO_tcr24" = WT_tcr24_vs_KO_tcr24[["up"]]
            )

pdf(paste0(out_dir_plots,"/02_WT_vs_KO_tcr24_naive_UPreg.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

s4 <- list( "down_WT_naive_vs_KO_naive" = WT_naive_vs_KO_naive[["down"]] ,
            "down_WT_tcr24_vs_KO_tcr24" = WT_tcr24_vs_KO_tcr24[["down"]]
            )

pdf(paste0(out_dir_plots,"/02_WT_vs_KO_tcr24_naive_DOWNreg.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()


##### multi-way:

s4 <- list( "up_WT_naive_to_TCR24" = WT_naive_to_TCR24[["up"]] ,
            "up_KO_naive_to_TCR24" = KO_naive_to_TCR24[["up"]],
            "up_WT_naive_vs_KO_naive" = WT_naive_vs_KO_naive[["up"]] ,
            "up_WT_tcr24_vs_KO_tcr24" = WT_tcr24_vs_KO_tcr24[["up"]],
            "down_WT_naive_vs_KO_naive" = WT_naive_vs_KO_naive[["down"]] ,
            "down_WT_tcr24_vs_KO_tcr24" = WT_tcr24_vs_KO_tcr24[["down"]]
            )

pdf(paste0(out_dir_plots,"/03_naive_to_TCR24_UP_vs_BL_CHANGES.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

s4 <- list( "down_WT_naive_to_TCR24" = WT_naive_to_TCR24[["down"]] ,
            "down_KO_naive_to_TCR24" = KO_naive_to_TCR24[["down"]],
            "up_WT_naive_vs_KO_naive" = WT_naive_vs_KO_naive[["up"]] ,
            "up_WT_tcr24_vs_KO_tcr24" = WT_tcr24_vs_KO_tcr24[["up"]],
            "down_WT_naive_vs_KO_naive" = WT_naive_vs_KO_naive[["down"]] ,
            "down_WT_tcr24_vs_KO_tcr24" = WT_tcr24_vs_KO_tcr24[["down"]]
            )

pdf(paste0(out_dir_plots,"/03_naive_to_TCR24_DOWN_vs_BL_CHANGES.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

ANY_change = unique(c(
    WT_naive_to_TCR24[["up"]],
    WT_naive_to_TCR24[["down"]],
    KO_naive_to_TCR24[["up"]] ,
    KO_naive_to_TCR24[["down"]],
    WT_naive_vs_KO_naive[["up"]], 
    WT_naive_vs_KO_naive[["down"]], 
    WT_tcr24_vs_KO_tcr24[["up"]] ,
    WT_tcr24_vs_KO_tcr24[["down"]]
)) 


mat_all = resALL[["norm_all"]][["Norm_expression"]]
mat_all = as.data.frame(mat_all[is.finite(rowSums(mat_all)),],stringsAsFactors=FALSE)
size = Ref_Anno_rn[rownames(mat_all),"width"] / 1000
mat_all = mat_all/size
mat_all[is.na(mat_all)] = 0
mat_all = mat_all[ANY_change ,c(   "WT_naive","WT_TCR24h", "KO_naive","KO_TCR24h")]
mat_all = as.data.frame(t(scale(t(mat_all),scale=TRUE,center=TRUE))[,colnames(mat_all)],stringsAsFactors=FALSE)

Plot = mat_all

sel_g = rep("ns",length( ANY_change ))
names(sel_g)= ANY_change

split_c = colnames(Plot)
split_c[grep("WT_naive",split_c)] = "01_WT"
split_c[grepl("WT_TCR24h",split_c)] = "01_WT"
split_c[grep("KO_naive",split_c)] = "02_KO"
split_c[grepl("KO_TCR24h",split_c)] = "02_KO"


sel_WT_naive_to_TCR24 = sel_g
sel_WT_naive_to_TCR24[names(sel_WT_naive_to_TCR24) %in%     WT_naive_to_TCR24[["up"]]] = "up"
sel_WT_naive_to_TCR24[names(sel_WT_naive_to_TCR24) %in%     WT_naive_to_TCR24[["down"]]] = "down"

sel_KO_naive_to_TCR24 = sel_g
sel_KO_naive_to_TCR24[names(sel_KO_naive_to_TCR24) %in%     KO_naive_to_TCR24[["up"]] ] = "up"
sel_KO_naive_to_TCR24[names(sel_KO_naive_to_TCR24) %in%     KO_naive_to_TCR24[["down"]]] = "down"
    
sel_WT_naive_vs_KO_naive = sel_g
sel_WT_naive_vs_KO_naive[names(sel_WT_naive_vs_KO_naive) %in%     WT_naive_vs_KO_naive[["up"]] ] = "up"
sel_WT_naive_vs_KO_naive[names(sel_WT_naive_vs_KO_naive) %in%     WT_naive_vs_KO_naive[["down"]]] = "down"

sel_WT_tcr24_vs_KO_tcr24 = sel_g
sel_WT_tcr24_vs_KO_tcr24[names(sel_WT_tcr24_vs_KO_tcr24) %in%     WT_tcr24_vs_KO_tcr24[["up"]] ] = "up"
sel_WT_tcr24_vs_KO_tcr24[names(sel_WT_tcr24_vs_KO_tcr24) %in%     WT_tcr24_vs_KO_tcr24[["down"]]] = "down"

ra= rowAnnotation(
                    sel_WT_naive_to_TCR24 = sel_WT_naive_to_TCR24 ,
                    sel_KO_naive_to_TCR24 = sel_KO_naive_to_TCR24 ,
                    sel_WT_naive_vs_KO_naive = sel_WT_naive_vs_KO_naive,
                    sel_WT_tcr24_vs_KO_tcr24 = sel_WT_tcr24_vs_KO_tcr24,

                    col = list( sel_WT_naive_to_TCR24 = c("ns" = "white", "up" = "red", "down" = "blue"),
                                sel_KO_naive_to_TCR24 = c("ns" = "white", "up" = "red", "down" = "blue"),
                                sel_WT_naive_vs_KO_naive = c("ns" = "white", "up" = "red", "down" = "blue"),
                                sel_WT_tcr24_vs_KO_tcr24 = c("ns" = "white", "up" = "red", "down" = "blue")  ),
                    
                    

                    Cyril_FCSinduced_sensitiveSRFlinked = anno_mark(at = which(rownames(Plot) %in% Cyril_FCSinduced_sensitiveSRFlinked), labels =rownames(Plot)[which(rownames(Plot) %in% Cyril_FCSinduced_sensitiveSRFlinked) ], labels_gp=gpar(fontsize = 1) ),
                    Gualdrini_SRFdirect = anno_mark(at = which(rownames(Plot) %in% Gualdrini_SRFdirect), labels =rownames(Plot)[which(rownames(Plot) %in% Gualdrini_SRFdirect) ], labels_gp=gpar(fontsize = 1) ),
                    annotation_name_gp = gpar(fontsize =3),
                    border=TRUE

                    
                )

mat_clust = Plot
mat_clust[is.na(mat_clust)] = 0


# MFA - UMAP:

GG = c( length(colnames(mat_all)[ grep("_naive",colnames(mat_all)) ]) ,length(colnames(mat_all)[ grep("_TCR24h",colnames(mat_all)) ]))
MFA <- MFA(     cbind(mat_clust[,grep("_naive",colnames(mat_all)) ],mat_clust[,grep("_TCR24h",colnames(mat_all)) ]),
                group = GG,type = rep("c",length(GG)),
                ncp=ncol(mat_clust),
                name.group=c("_naive","_TCR24h"),
                graph=FALSE)

indiv = MFA$ind$coord

# rpca = FactoMineR::PCA(mat_clust, graph = FALSE,scale.unit = FALSE,ncp = ncol(mat_clust))
# indiv = rpca$ind$coord

k = round(( dim(indiv)[1]   )^(1/2)  ) 

knn.real = get.knn(as.matrix(indiv), k = k, algorithm="brute")
knn.real = data.frame(from = rep(1:nrow(knn.real$nn.index), k), to = as.vector(knn.real$nn.index), weight = 1/(1 + as.vector(knn.real$nn.dist)))
nw.real = graph_from_data_frame(knn.real, directed = FALSE)
nw.real = simplify(nw.real)


#    cluster_leading_eigen cluster_louvain cluster_walktrap
lc.real = cluster_louvain(nw.real) 
louvain = as.factor(membership(lc.real)) 
cl = paste0("Cluster_",louvain)
names(cl) = rownames(indiv)
cln=cl
HT = Heatmap(     Plot,
                        width = unit(20, "mm"),
                        row_split = paste0("Cluster_",cl) , 
                        column_split = split_c,
                        column_gap = unit(1, "mm"),
                        row_gap = unit(1, "mm"),
                        cluster_column_slices = FALSE,
                        cluster_columns=FALSE, 
                        right_annotation = ra,
                        show_row_names = FALSE,
                        show_column_names = TRUE,
                        column_title_gp = gpar(fontsize = 3),
                        row_title_gp = gpar(fontsize = 3),
                        heatmap_legend_param = list(title = "z-scored ATACsignal",direction = "horizontal"),
                        column_names_gp = gpar(fontsize = 2),
                        row_dend_width = unit(2, "mm"),
                        show_row_dend = TRUE,
                        row_dend_reorder = TRUE,
                        column_dend_reorder = FALSE,
                        use_raster = TRUE,
                        na_col = "white",
                        border=TRUE,
                        col = colorRamp2( seq(-2, 2, by = 0.1), cividis(length(seq(-2, 2, by = 0.1))))
                        )
pdf(paste0(out_dir_plots,"/04_Heat_all_intron.pdf"),useDingbats = FALSE)
    draw(HT, heatmap_legend_side="right")
dev.off()

Plot$Clusters = NA
Plot[names(cln) , "Clusters" ] = cln
write.table(Plot,paste0(out_dir_plots,"/04_Heat_Clusters.txt"),sep="\t",col.names=NA)
Plot_04 = Plot


cl = cln[names(cln) %in% Universe]

list_comp = split(names(cl),cl)
check = unique(unlist(list_comp))
Universe = union(Universe,check)

msig_plots = paste0(out_dir_plots,"/msig_HJ_NAIVE_TCR24h/")
unlink(msig_plots, recursive=TRUE)
dir.create(msig_plots)

# for each list assess the overlap:
# HG + JI:
HGJa = lapply( 1:length(MSIG_list_symbol_entrez_mouse), function(j)  {

     l1 = lapply(1:length(list_comp),function(lv) {
        l2 = lapply( 1:length(MSIG_list_symbol_entrez_mouse[[j]]), function(y)  {
            set1=list_comp[[lv]] #gsub(".*_","",list_comp[[lv]])
            set2=MSIG_list_symbol_entrez_mouse[[j]][[y]]#gsub(".*_","",MSIG_list_symbol_entrez_mouse[[j]][[y]])

            HGu=HG_fun_upper(set1,set2,Universe)
            HGl=HG_fun_lower(set1,set2,Universe)
            HG = min(c(HGu,HGl))

            JI= JI_fun(set1,set2)

            Genes = paste0(intersect(set1,set2),collapse=";")

            j_nm<-c()
            for (nm_rep in 1:500){
                A_<-sample(Universe,length(set1))
                p_rand = JI_fun(A_,set2)
                j_nm<-c(j_nm,p_rand)
            }

            SCALED_JI<-(JI-mean(j_nm))/sd(j_nm)
            if(is.nan(SCALED_JI)){SCALED_JI=0}
            return(c( HG=HG,HGu=HGu,HGl=HGl,JI= JI, SCALED_JI = SCALED_JI,Genes=Genes ))

        })
        names(l2) = names(MSIG_list_symbol_entrez_mouse[[j]])
        l2 = as.data.frame(do.call(rbind,l2),stringsAsFactors=FALSE)
        HG = as.numeric(as.character(l2$HG)) #
        names(HG) = rownames(l2)        
        HGu = as.numeric(as.character(l2$HGu)) #
        names(HGu) = rownames(l2)
        HGl = as.numeric(as.character(l2$HGl)) #
        names(HGl) = rownames(l2)
        JI = as.numeric(as.character(l2$JI))
        names(JI) = rownames(l2)
        SCALED_JI = as.numeric(as.character(l2$SCALED_JI))
        names(SCALED_JI) = rownames(l2)
        Genes = l2$Genes
        names(Genes) = rownames(l2)
        return(list(HG =HG, HGu = HGu,HGl=HGl,JI=JI,SCALED_JI=SCALED_JI,Genes=Genes))
    })

    names(l1) = names(list_comp)
    HG = as.data.frame(do.call(cbind,lapply(l1, `[[`, 1)),stringsAsFactors=FALSE)    
    HGu = as.data.frame(do.call(cbind,lapply(l1, `[[`, 2)),stringsAsFactors=FALSE)
    HGl = as.data.frame(do.call(cbind,lapply(l1, `[[`, 3)),stringsAsFactors=FALSE)
    JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 4)),stringsAsFactors=FALSE)
    SCALED_JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 5)),stringsAsFactors=FALSE)
    Genes = as.data.frame(do.call(cbind,lapply(l1, `[[`, 6)),stringsAsFactors=FALSE)

    res_l = list(   
                    HG = as.data.frame(do.call(cbind,lapply(l1, `[[`, 1)),stringsAsFactors=FALSE)    ,
                    HGu = as.data.frame(do.call(cbind,lapply(l1, `[[`, 2)),stringsAsFactors=FALSE),
                    HGl = as.data.frame(do.call(cbind,lapply(l1, `[[`, 3)),stringsAsFactors=FALSE),
                    JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 4)),stringsAsFactors=FALSE),
                    SCALED_JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 5)),stringsAsFactors=FALSE),
                    FDRu = apply(HGu,2,function(x){p.adjust(x,method='BH')}),
                    FDRl = apply(HGl,2,function(x){p.adjust(x,method='BH')}),
                    Genes = as.data.frame(do.call(cbind,lapply(l1, `[[`, 6)),stringsAsFactors=FALSE)
                    )

    output_file <- paste0(msig_plots,"/",names(MSIG_list_symbol_entrez_mouse)[j], "_MSIG_Enrichment_Results.xlsx")
    write.xlsx( res_l, output_file ,row.names = TRUE )

    return(list(HG=HG,HGu = HGu,HGl=HGl,JI=JI,SCALED_JI=SCALED_JI,Genes=Genes))

})

names(HGJa) = names(MSIG_list_symbol_entrez_mouse)


for(xx in names(MSIG_list_symbol_entrez_mouse)){

    HGT = HGJa[[xx]][["HG"]]
    HGT = apply(HGT,2,function(x){p.adjust(x,method='BH')})
    HGT = as.data.frame(HGT,stringsAsFactors=FALSE)
    HGT=-log10(HGT)
    max_replace = max(c(max(unlist(HGT)[is.finite(unlist(HGT))])))
    HGT[] <- lapply(HGT, function(i) if(is.numeric(i)) ifelse(is.infinite(i), max_replace, i) else i)
    ww = which(apply(HGT, MARGIN=c(1), max) >= 2   )
    HGT = HGT[ww,,drop=FALSE]
    compress = 3
    HGT[HGT>=compress] = compress
    HGTrad = sqrt(HGT/pi)
    HGTrad = HGTrad/max(HGTrad)

    SCALED_JI = HGJa[[xx]][["SCALED_JI"]]


    rr = intersect(rownames(SCALED_JI),rownames(HGTrad))
    cc = intersect(colnames(SCALED_JI),colnames(HGTrad))

    SCALED_JI = SCALED_JI[rr,cc,drop=FALSE]
    HGTrad = HGTrad[rr,cc,drop=FALSE]

    if(dim(HGTrad)[1] !=0 & dim(HGTrad)[2] !=0){

        col_fun = colorRamp2(c(-1*max(abs(SCALED_JI)),-1*max(abs(SCALED_JI))/2,0,max(abs(SCALED_JI))/2,max(abs(SCALED_JI))),c("#034961","#51D06D", "white","#FF9E39","#C51300")  )
        
        htord = Heatmap(    SCALED_JI, 
                                cluster_rows = TRUE, 
                                row_dend_reorder = TRUE,
                                cluster_columns = TRUE,
                                column_dend_reorder=TRUE
                            )
        cO = colnames(SCALED_JI)[unlist(column_order(htord))]
        rO = rownames(SCALED_JI)[unlist(row_order(htord))]
        SCALED_JI = SCALED_JI[rO,cO,drop=FALSE]
        HGTrad = HGTrad[rownames(SCALED_JI),colnames(SCALED_JI),drop=FALSE]
        dds = 20/max(nrow(SCALED_JI),ncol(SCALED_JI))
        
        HR1 = Heatmap(          SCALED_JI,

                                column_title = paste0("Fraction: radius min: ",round(min(HGJa[[xx]][["HG"]]),4)," radius max: ",round(max(HGJa[[xx]][["HG"]]),4)),
                                width = ncol(SCALED_JI)*unit(dds, "cm"), 
                                height = nrow(SCALED_JI)*unit(dds, "cm"),
                                #rect_gp = gpar(type = "none"), 
                                #cell_fun = function(j, i, x, y, width, height, fill) {
                                #    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = NA, fill = NA))
                                #    grid.circle(x = x, y = y, r = abs(HGTrad[i, j])/2.5 * unit(dds, "cm"), 
                                #    gp = gpar(fill = col_fun(SCALED_JI[i, j]), col = NA))
                                #}, 
                                show_row_names = TRUE,
                                show_column_names = TRUE,
                                column_title_gp = gpar(fontsize = 2),
                                row_title_gp = gpar(fontsize = 1),
                                row_names_gp = gpar(fontsize = 1),
                                heatmap_legend_param = list(title = "Calls",direction = "horizontal"),
                                column_names_gp = gpar(fontsize = 2),
                                row_dend_width = unit(2, "mm"),
                                border=TRUE,
                                col = col_fun ,
                                cluster_rows = TRUE, 
                                cluster_columns = TRUE,
                                row_dend_reorder = TRUE,
                                column_dend_reorder=TRUE

                                )
        pdf(paste0(msig_plots,"/", gsub(" |:","_",xx) ,"_enrichment.pdf"),useDingbats = FALSE , width=10   , height=10)
                draw(HR1)
        dev.off()  
    }
}

####  TCR24 RestedIL12 #########
# Extract the differences ######
################################
# extract changes that are transition specific:

WT_restedIL12_vs_WT_TCR24h = list()
up_all_WT_restedIL12_vs_WT_TCR24h <- rownames(resALL[["norm_all"]][["WT_restedIL12_vs_WT_TCR24h"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_all_WT_restedIL12_vs_WT_TCR24h <-  rownames(resALL[["norm_all"]][["WT_restedIL12_vs_WT_TCR24h"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_WT_restedIL12_vs_WT_TCR24h <- rownames(resALL[["norm_intron"]][["WT_restedIL12_vs_WT_TCR24h"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_WT_restedIL12_vs_WT_TCR24h <-  rownames(resALL[["norm_intron"]][["WT_restedIL12_vs_WT_TCR24h"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
WT_restedIL12_vs_WT_TCR24h[["up"]] = union(up_all_WT_restedIL12_vs_WT_TCR24h,up_intron_WT_restedIL12_vs_WT_TCR24h)
WT_restedIL12_vs_WT_TCR24h[["down"]] = union(down_all_WT_restedIL12_vs_WT_TCR24h,down_intron_WT_restedIL12_vs_WT_TCR24h)

KO_restedIL12_vs_KO_TCR24h = list()
up_all_KO_restedIL12_vs_KO_TCR24h <- rownames(resALL[["norm_all"]][["KO_restedIL12_vs_KO_TCR24h"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_all_KO_restedIL12_vs_KO_TCR24h <-  rownames(resALL[["norm_all"]][["KO_restedIL12_vs_KO_TCR24h"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_KO_restedIL12_vs_KO_TCR24h <- rownames(resALL[["norm_intron"]][["KO_restedIL12_vs_KO_TCR24h"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_KO_restedIL12_vs_KO_TCR24h <-  rownames(resALL[["norm_intron"]][["KO_restedIL12_vs_KO_TCR24h"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
KO_restedIL12_vs_KO_TCR24h[["up"]] = union(up_all_KO_restedIL12_vs_KO_TCR24h,up_intron_KO_restedIL12_vs_KO_TCR24h)
KO_restedIL12_vs_KO_TCR24h[["down"]] = union(down_all_KO_restedIL12_vs_KO_TCR24h,down_intron_KO_restedIL12_vs_KO_TCR24h)

s4 <- list( "up_WT_restedIL12_vs_WT_TCR24h" = WT_restedIL12_vs_WT_TCR24h[["up"]] ,
            "up_KO_restedIL12_vs_KO_TCR24h" = KO_restedIL12_vs_KO_TCR24h[["up"]]
            )

pdf(paste0(out_dir_plots,"/05_tcr_rested_transition_up_in_WT_KO.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

s4 <- list( "down_WT_restedIL12_vs_WT_TCR24h" = WT_restedIL12_vs_WT_TCR24h[["down"]] ,
            "down_KO_restedIL12_vs_KO_TCR24h" = KO_restedIL12_vs_KO_TCR24h[["down"]]
            )

pdf(paste0(out_dir_plots,"/05_tcr_rested_transition_down_in_WT_KO.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

KO_restedIL12_vs_WT_restedIL12 = list()
up_all_KO_restedIL12_vs_WT_restedIL12 <- rownames(resALL[["norm_all"]][["KO_restedIL12_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_all_KO_restedIL12_vs_WT_restedIL12 <-  rownames(resALL[["norm_all"]][["KO_restedIL12_vs_WT_restedIL12"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_KO_restedIL12_vs_WT_restedIL12 <- rownames(resALL[["norm_intron"]][["KO_restedIL12_vs_WT_restedIL12"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_KO_restedIL12_vs_WT_restedIL12 <-  rownames(resALL[["norm_intron"]][["KO_restedIL12_vs_WT_restedIL12"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
KO_restedIL12_vs_WT_restedIL12[["up"]] = union(up_all_KO_restedIL12_vs_WT_restedIL12,up_intron_KO_restedIL12_vs_WT_restedIL12)
KO_restedIL12_vs_WT_restedIL12[["down"]] = union(down_all_KO_restedIL12_vs_WT_restedIL12,down_intron_KO_restedIL12_vs_WT_restedIL12)

KO_TCR24h_vs_WT_TCR24h = list()
up_all_KO_TCR24h_vs_WT_TCR24h <- rownames(resALL[["norm_all"]][["KO_TCR24h_vs_WT_TCR24h"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_all_KO_TCR24h_vs_WT_TCR24h <-  rownames(resALL[["norm_all"]][["KO_TCR24h_vs_WT_TCR24h"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_KO_TCR24h_vs_WT_TCR24h <- rownames(resALL[["norm_intron"]][["KO_TCR24h_vs_WT_TCR24h"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_KO_TCR24h_vs_WT_TCR24h <-  rownames(resALL[["norm_intron"]][["KO_TCR24h_vs_WT_TCR24h"]] %>% dplyr::filter(as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
KO_TCR24h_vs_WT_TCR24h[["up"]] = union(up_all_KO_TCR24h_vs_WT_TCR24h,up_intron_KO_TCR24h_vs_WT_TCR24h)
KO_TCR24h_vs_WT_TCR24h[["down"]] = union(down_all_KO_TCR24h_vs_WT_TCR24h,down_intron_KO_TCR24h_vs_WT_TCR24h)

s4 <- list( "up_KO_restedIL12_vs_WT_restedIL12" = KO_restedIL12_vs_WT_restedIL12[["up"]] ,
            "up_KO_TCR24h_vs_WT_TCR24h" = KO_TCR24h_vs_WT_TCR24h[["up"]]
            )

pdf(paste0(out_dir_plots,"/06_WT_vs_KO_tcr24_rested_UPreg.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

s4 <- list( "down_KO_restedIL12_vs_WT_restedIL12" = KO_restedIL12_vs_WT_restedIL12[["down"]] ,
            "down_KO_TCR24h_vs_WT_TCR24h" = KO_TCR24h_vs_WT_TCR24h[["down"]]
            )

pdf(paste0(out_dir_plots,"/06_WT_vs_KO_tcr24_rested_DOWNreg.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()


##### multi-way:

s4 <- list( "up_WT_restedIL12_vs_WT_TCR24h" = WT_restedIL12_vs_WT_TCR24h[["up"]] ,
            "up_KO_restedIL12_vs_KO_TCR24h" = KO_restedIL12_vs_KO_TCR24h[["up"]],
            "up_KO_restedIL12_vs_WT_restedIL12" = KO_restedIL12_vs_WT_restedIL12[["up"]] ,
            "up_KO_TCR24h_vs_WT_TCR24h" = KO_TCR24h_vs_WT_TCR24h[["up"]],
            "down_KO_restedIL12_vs_WT_restedIL12" = KO_restedIL12_vs_WT_restedIL12[["down"]] ,
            "down_KO_TCR24h_vs_WT_TCR24h" = KO_TCR24h_vs_WT_TCR24h[["down"]]
            )

pdf(paste0(out_dir_plots,"/07_TCR24_to_redted_UP_vs_BL_CHANGES.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

s4 <- list( "down_WT_restedIL12_vs_WT_TCR24h" = WT_restedIL12_vs_WT_TCR24h[["down"]] ,
            "down_KO_restedIL12_vs_KO_TCR24h" = KO_restedIL12_vs_KO_TCR24h[["down"]],
            "up_KO_restedIL12_vs_WT_restedIL12" = KO_restedIL12_vs_WT_restedIL12[["up"]] ,
            "up_KO_TCR24h_vs_WT_TCR24h" = KO_TCR24h_vs_WT_TCR24h[["up"]],
            "down_KO_restedIL12_vs_WT_restedIL12" = KO_restedIL12_vs_WT_restedIL12[["down"]] ,
            "down_KO_TCR24h_vs_WT_TCR24h" = KO_TCR24h_vs_WT_TCR24h[["down"]]
            )

pdf(paste0(out_dir_plots,"/07_TCR24_to_redted_DOWN_vs_BL_CHANGES.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

ANY_change = unique(c(
    WT_restedIL12_vs_WT_TCR24h[["up"]],
    WT_restedIL12_vs_WT_TCR24h[["down"]],
    KO_restedIL12_vs_KO_TCR24h[["up"]] ,
    KO_restedIL12_vs_KO_TCR24h[["down"]],
    KO_restedIL12_vs_WT_restedIL12[["up"]], 
    KO_restedIL12_vs_WT_restedIL12[["down"]], 
    KO_TCR24h_vs_WT_TCR24h[["up"]] ,
    KO_TCR24h_vs_WT_TCR24h[["down"]]
)) 


mat_all = resALL[["norm_all"]][["Norm_expression"]]
mat_all = as.data.frame(mat_all[is.finite(rowSums(mat_all)),],stringsAsFactors=FALSE)
size = Ref_Anno_rn[rownames(mat_all),"width"] / 1000
mat_all = mat_all/size
mat_all[is.na(mat_all)] = 0
mat_all = mat_all[ANY_change ,c( "WT_TCR24h","WT_restedIL12","KO_TCR24h","KO_restedIL12")]
mat_all = as.data.frame(t(scale(t(mat_all),scale=TRUE,center=TRUE))[,colnames(mat_all)],stringsAsFactors=FALSE)

Plot = mat_all

sel_g = rep("ns",length( ANY_change ))
names(sel_g)= ANY_change

split_c = colnames(Plot)
split_c[grepl("WT_TCR24h",split_c)] = "01_WT"
split_c[grep("WT_restedIL12",split_c)] = "01_WT"
split_c[grepl("KO_TCR24h",split_c)] = "02_KO"
split_c[grep("KO_restedIL12",split_c)] = "02_KO"

sel_WT_restedIL12_vs_WT_TCR24h = sel_g
sel_WT_restedIL12_vs_WT_TCR24h[names(sel_WT_restedIL12_vs_WT_TCR24h) %in%     WT_restedIL12_vs_WT_TCR24h[["up"]]] = "up"
sel_WT_restedIL12_vs_WT_TCR24h[names(sel_WT_restedIL12_vs_WT_TCR24h) %in%     WT_restedIL12_vs_WT_TCR24h[["down"]]] = "down"

sel_KO_restedIL12_vs_KO_TCR24h = sel_g
sel_KO_restedIL12_vs_KO_TCR24h[names(sel_KO_restedIL12_vs_KO_TCR24h) %in%     KO_restedIL12_vs_KO_TCR24h[["up"]] ] = "up"
sel_KO_restedIL12_vs_KO_TCR24h[names(sel_KO_restedIL12_vs_KO_TCR24h) %in%     KO_restedIL12_vs_KO_TCR24h[["down"]]] = "down"
    
sel_KO_restedIL12_vs_WT_restedIL12 = sel_g
sel_KO_restedIL12_vs_WT_restedIL12[names(sel_KO_restedIL12_vs_WT_restedIL12) %in%     KO_restedIL12_vs_WT_restedIL12[["up"]] ] = "up"
sel_KO_restedIL12_vs_WT_restedIL12[names(sel_KO_restedIL12_vs_WT_restedIL12) %in%     KO_restedIL12_vs_WT_restedIL12[["down"]]] = "down"

sel_KO_TCR24h_vs_WT_TCR24h = sel_g
sel_KO_TCR24h_vs_WT_TCR24h[names(sel_KO_TCR24h_vs_WT_TCR24h) %in%     KO_TCR24h_vs_WT_TCR24h[["up"]] ] = "up"
sel_KO_TCR24h_vs_WT_TCR24h[names(sel_KO_TCR24h_vs_WT_TCR24h) %in%     KO_TCR24h_vs_WT_TCR24h[["down"]]] = "down"

ra= rowAnnotation(
                    sel_WT_restedIL12_vs_WT_TCR24h = sel_WT_restedIL12_vs_WT_TCR24h ,
                    sel_KO_restedIL12_vs_KO_TCR24h = sel_KO_restedIL12_vs_KO_TCR24h ,
                    sel_KO_restedIL12_vs_WT_restedIL12 = sel_KO_restedIL12_vs_WT_restedIL12,
                    sel_KO_TCR24h_vs_WT_TCR24h = sel_KO_TCR24h_vs_WT_TCR24h,

                    col = list( sel_WT_restedIL12_vs_WT_TCR24h = c("ns" = "white", "up" = "red", "down" = "blue"),
                                sel_KO_restedIL12_vs_KO_TCR24h = c("ns" = "white", "up" = "red", "down" = "blue"),
                                sel_KO_restedIL12_vs_WT_restedIL12 = c("ns" = "white", "up" = "red", "down" = "blue"),
                                sel_KO_TCR24h_vs_WT_TCR24h = c("ns" = "white", "up" = "red", "down" = "blue")  ),                    
                    Cyril_FCSinduced_sensitiveSRFlinked = anno_mark(at = which(rownames(Plot) %in% Cyril_FCSinduced_sensitiveSRFlinked), labels =rownames(Plot)[which(rownames(Plot) %in% Cyril_FCSinduced_sensitiveSRFlinked) ], labels_gp=gpar(fontsize = 1) ),
                    Gualdrini_SRFdirect = anno_mark(at = which(rownames(Plot) %in% Gualdrini_SRFdirect), labels =rownames(Plot)[which(rownames(Plot) %in% Gualdrini_SRFdirect) ], labels_gp=gpar(fontsize = 1) ),
                    annotation_name_gp = gpar(fontsize =3),
                    border=TRUE
                )

mat_clust = Plot
mat_clust[is.na(mat_clust)] = 0

# MFA - UMAP:

GG = c( length(colnames(mat_clust)[ grep("_TCR24h",colnames(mat_clust)) ]) ,length(colnames(mat_clust)[ grep("_restedIL12",colnames(mat_clust)) ]))
MFA <- MFA(     cbind(mat_clust[,grep("_TCR24h",colnames(mat_clust)) ],mat_clust[,grep("_restedIL12",colnames(mat_clust)) ]),
                group = GG,type = rep("c",length(GG)),
                ncp=ncol(mat_clust),
                name.group=c("_TCR24h","_restedIL12"),
                graph=FALSE)
indiv = MFA$ind$coord

# rpca = FactoMineR::PCA(mat_clust, graph = FALSE,scale.unit = FALSE,ncp = ncol(mat_clust))
# indiv = rpca$ind$coord

k = round(( dim(indiv)[1]   )^(1/2) ) 

knn.real = get.knn(as.matrix(indiv), k = k, algorithm="brute")
knn.real = data.frame(from = rep(1:nrow(knn.real$nn.index), k), to = as.vector(knn.real$nn.index), weight = 1/(1 + as.vector(knn.real$nn.dist)))
nw.real = graph_from_data_frame(knn.real, directed = FALSE)
nw.real = simplify(nw.real)
lc.real = cluster_louvain(nw.real)
louvain = as.factor(membership(lc.real)) 

# cl = kmeans(mat_clust, centers = 10)$cluster

cl = paste0("Cluster_",louvain)
names(cl) = rownames(indiv)
cln=cl

HT = Heatmap(     Plot,
                        width = unit(20, "mm"),
                        row_split = paste0("Cluster_",cl) , 
                        column_split = split_c,
                        column_gap = unit(1, "mm"),
                        row_gap = unit(1, "mm"),
                        cluster_column_slices = FALSE,
                        cluster_columns=FALSE, 
                        right_annotation = ra,
                        show_row_names = FALSE,
                        show_column_names = TRUE,
                        column_title_gp = gpar(fontsize = 3),
                        row_title_gp = gpar(fontsize = 3),
                        heatmap_legend_param = list(title = "z-scored ATACsignal",direction = "horizontal"),
                        column_names_gp = gpar(fontsize = 2),
                        row_dend_width = unit(2, "mm"),
                        show_row_dend = TRUE,
                                            use_raster = TRUE,
                        row_dend_reorder = TRUE,
                        column_dend_reorder = FALSE,
                        na_col = "white",
                        border=TRUE,
                        col = colorRamp2( seq(-2, 2, by = 0.1), cividis(length(seq(-2, 2, by = 0.1))))
                        )
pdf(paste0(out_dir_plots,"/08_Heat_all_intron.pdf"),useDingbats = FALSE)
    draw(HT, heatmap_legend_side="right")
dev.off()

Plot$Clusters = NA
Plot[names(cln) , "Clusters" ] = cln
write.table(Plot,paste0(out_dir_plots,"/08_Heat_Clusters.txt"),sep="\t",col.names=NA)
Plot_08 = Plot

cl = cln[names(cln) %in% Universe]
list_comp = split(names(cl),cl)
check = unique(unlist(list_comp))
Universe = union(Universe,check)

# for each list assess the overlap:
# HG + JI:

msig_plots = paste0(out_dir_plots,"/msig_HJ_TCR24h_RESTED/")
unlink(msig_plots, recursive=TRUE)
dir.create(msig_plots)

# for each list assess the overlap:
# HG + JI:
HGJb = lapply( 1:length(MSIG_list_symbol_entrez_mouse), function(j)  {

     l1 = lapply(1:length(list_comp),function(lv) {
        l2 = lapply( 1:length(MSIG_list_symbol_entrez_mouse[[j]]), function(y)  {
            set1=list_comp[[lv]] #gsub(".*_","",list_comp[[lv]])
            set2=MSIG_list_symbol_entrez_mouse[[j]][[y]]#gsub(".*_","",MSIG_list_symbol_entrez_mouse[[j]][[y]])

            HGu=HG_fun_upper(set1,set2,Universe)
            HGl=HG_fun_lower(set1,set2,Universe)
            HG = min(c(HGu,HGl))

            JI= JI_fun(set1,set2)

            Genes = paste0(intersect(set1,set2),collapse=";")

            j_nm<-c()
            for (nm_rep in 1:500){
                A_<-sample(Universe,length(set1))
                p_rand = JI_fun(A_,set2)
                j_nm<-c(j_nm,p_rand)
            }

            SCALED_JI<-(JI-mean(j_nm))/sd(j_nm)
            if(is.nan(SCALED_JI)){SCALED_JI=0}
            return(c( HG=HG,HGu=HGu,HGl=HGl,JI= JI, SCALED_JI = SCALED_JI,Genes=Genes ))

        })
        names(l2) = names(MSIG_list_symbol_entrez_mouse[[j]])
        l2 = as.data.frame(do.call(rbind,l2),stringsAsFactors=FALSE)
        HG = as.numeric(as.character(l2$HG)) #
        names(HG) = rownames(l2)        
        HGu = as.numeric(as.character(l2$HGu)) #
        names(HGu) = rownames(l2)
        HGl = as.numeric(as.character(l2$HGl)) #
        names(HGl) = rownames(l2)
        JI = as.numeric(as.character(l2$JI))
        names(JI) = rownames(l2)
        SCALED_JI = as.numeric(as.character(l2$SCALED_JI))
        names(SCALED_JI) = rownames(l2)
        Genes = l2$Genes
        names(Genes) = rownames(l2)
        return(list(HG =HG, HGu = HGu,HGl=HGl,JI=JI,SCALED_JI=SCALED_JI,Genes=Genes))
    })

    names(l1) = names(list_comp)
    HG = as.data.frame(do.call(cbind,lapply(l1, `[[`, 1)),stringsAsFactors=FALSE)    
    HGu = as.data.frame(do.call(cbind,lapply(l1, `[[`, 2)),stringsAsFactors=FALSE)
    HGl = as.data.frame(do.call(cbind,lapply(l1, `[[`, 3)),stringsAsFactors=FALSE)
    JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 4)),stringsAsFactors=FALSE)
    SCALED_JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 5)),stringsAsFactors=FALSE)
    Genes = as.data.frame(do.call(cbind,lapply(l1, `[[`, 6)),stringsAsFactors=FALSE)

    res_l = list(   
                    HG = as.data.frame(do.call(cbind,lapply(l1, `[[`, 1)),stringsAsFactors=FALSE)    ,
                    HGu = as.data.frame(do.call(cbind,lapply(l1, `[[`, 2)),stringsAsFactors=FALSE),
                    HGl = as.data.frame(do.call(cbind,lapply(l1, `[[`, 3)),stringsAsFactors=FALSE),
                    JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 4)),stringsAsFactors=FALSE),
                    SCALED_JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 5)),stringsAsFactors=FALSE),
                    FDRu = apply(HGu,2,function(x){p.adjust(x,method='BH')}),
                    FDRl = apply(HGl,2,function(x){p.adjust(x,method='BH')}),
                    Genes = as.data.frame(do.call(cbind,lapply(l1, `[[`, 6)),stringsAsFactors=FALSE)
                    )

    output_file <- paste0(msig_plots,"/",names(MSIG_list_symbol_entrez_mouse)[j], "_MSIG_Enrichment_Results.xlsx")
    write.xlsx( res_l, output_file ,row.names = TRUE )

    return(list(HG=HG,HGu = HGu,HGl=HGl,JI=JI,SCALED_JI=SCALED_JI,Genes=Genes))

})


names(HGJb) = names(MSIG_list_symbol_entrez_mouse)


for(xx in names(MSIG_list_symbol_entrez_mouse)){

    HGT = HGJb[[xx]][["HG"]]
    HGT = apply(HGT,2,function(x){p.adjust(x,method='BH')})
    HGT = as.data.frame(HGT,stringsAsFactors=FALSE)
    HGT=-log10(HGT)
    max_replace = max(c(max(unlist(HGT)[is.finite(unlist(HGT))])))
    HGT[] <- lapply(HGT, function(i) if(is.numeric(i)) ifelse(is.infinite(i), max_replace, i) else i)
    ww = which(apply(HGT, MARGIN=c(1), max) >= 2   )
    HGT = HGT[ww,,drop=FALSE]
    compress = 3
    HGT[HGT>=compress] = compress
    HGTrad = sqrt(HGT/pi)
    HGTrad = HGTrad/max(HGTrad)

    SCALED_JI = HGJb[[xx]][["SCALED_JI"]]

    rr = intersect(rownames(SCALED_JI),rownames(HGTrad))
    cc = intersect(colnames(SCALED_JI),colnames(HGTrad))

    SCALED_JI = SCALED_JI[rr,cc,drop=FALSE]
    HGTrad = HGTrad[rr,cc,drop=FALSE]

    if(dim(HGTrad)[1] !=0 & dim(HGTrad)[2] !=0){

        col_fun = colorRamp2(c(-1*max(abs(SCALED_JI)),-1*max(abs(SCALED_JI))/2,0,max(abs(SCALED_JI))/2,max(abs(SCALED_JI))),c("#034961","#51D06D", "white","#FF9E39","#C51300")  )
        
        htord = Heatmap(    SCALED_JI, 
                                cluster_rows = TRUE, 
                                row_dend_reorder = TRUE,
                                cluster_columns = TRUE,
                                column_dend_reorder=TRUE
                            )
        cO = colnames(SCALED_JI)[unlist(column_order(htord))]
        rO = rownames(SCALED_JI)[unlist(row_order(htord))]
        SCALED_JI = SCALED_JI[rO,cO,drop=FALSE]
        HGTrad = HGTrad[rownames(SCALED_JI),colnames(SCALED_JI),drop=FALSE]
        dds = 20/max(nrow(SCALED_JI),ncol(SCALED_JI))
               HR1 = Heatmap(          SCALED_JI,

                                column_title = paste0("Fraction: radius min: ",round(min(HGJb[[xx]][["HG"]]),4)," radius max: ",round(max(HGJb[[xx]][["HG"]]),4)),
                                width = ncol(SCALED_JI)*unit(dds, "cm"), 
                                height = nrow(SCALED_JI)*unit(dds, "cm"),
                                #rect_gp = gpar(type = "none"), 
                                #cell_fun = function(j, i, x, y, width, height, fill) {
                                #    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = NA, fill = NA))
                                #    grid.circle(x = x, y = y, r = abs(HGTrad[i, j])/2.5 * unit(dds, "cm"), 
                                #    gp = gpar(fill = col_fun(SCALED_JI[i, j]), col = NA))
                                #}, 
                                show_row_names = TRUE,
                                show_column_names = TRUE,
                                column_title_gp = gpar(fontsize = 2),
                                row_title_gp = gpar(fontsize = 1),
                                row_names_gp = gpar(fontsize = 1),
                                heatmap_legend_param = list(title = "Calls",direction = "horizontal"),
                                column_names_gp = gpar(fontsize = 2),
                                row_dend_width = unit(2, "mm"),
                                border=TRUE,
                                col = col_fun ,
                                cluster_rows = TRUE, 
                                cluster_columns = TRUE,
                                row_dend_reorder = TRUE,
                                column_dend_reorder=TRUE

                                )
        pdf(paste0(msig_plots,"/", gsub(" |:","_",xx) ,"_enrichment.pdf"),useDingbats = FALSE , width=10   , height=10)
                draw(HR1)
        dev.off()  
    }
}


################################################################################################################################
################################################################################################################################

#### IL2 TC analysis############
# Extract the differences ######
################################
# extract changes that are transition specific:


WT_restedIL12_vs_individual_TP = list()
# up
up_all_WT_30minIL2_vs_WT_restedIL12 <- rownames(resALL[["norm_all"]][["WT_30minIL2_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_all_WT_60minIL2_vs_WT_restedIL12 <- rownames(resALL[["norm_all"]][["WT_60minIL2_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_all_WT_120minIL2_vs_WT_restedIL12 <- rownames(resALL[["norm_all"]][["WT_120minIL2_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
# down
down_all_WT_30minIL2_vs_WT_restedIL12 <- rownames(resALL[["norm_all"]][["WT_30minIL2_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_all_WT_60minIL2_vs_WT_restedIL12 <- rownames(resALL[["norm_all"]][["WT_60minIL2_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_all_WT_120minIL2_vs_WT_restedIL12 <- rownames(resALL[["norm_all"]][["WT_120minIL2_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))

# up
up_intron_WT_30minIL2_vs_WT_restedIL12 <- rownames(resALL[["norm_intron"]][["WT_30minIL2_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_WT_60minIL2_vs_WT_restedIL12 <- rownames(resALL[["norm_intron"]][["WT_60minIL2_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_WT_120minIL2_vs_WT_restedIL12 <- rownames(resALL[["norm_intron"]][["WT_120minIL2_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
# down
down_intron_WT_30minIL2_vs_WT_restedIL12 <- rownames(resALL[["norm_intron"]][["WT_30minIL2_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_WT_60minIL2_vs_WT_restedIL12 <- rownames(resALL[["norm_intron"]][["WT_60minIL2_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_WT_120minIL2_vs_WT_restedIL12 <- rownames(resALL[["norm_intron"]][["WT_120minIL2_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))

WT_restedIL12_vs_individual_TP[["up_30"]] = union(up_all_WT_30minIL2_vs_WT_restedIL12,up_intron_WT_30minIL2_vs_WT_restedIL12)
WT_restedIL12_vs_individual_TP[["up_60"]] = union(up_all_WT_60minIL2_vs_WT_restedIL12,up_intron_WT_60minIL2_vs_WT_restedIL12)
WT_restedIL12_vs_individual_TP[["up_120"]] = union(up_all_WT_120minIL2_vs_WT_restedIL12,up_intron_WT_120minIL2_vs_WT_restedIL12)

WT_restedIL12_vs_individual_TP[["down_30"]] = union(down_all_WT_30minIL2_vs_WT_restedIL12,down_intron_WT_30minIL2_vs_WT_restedIL12)
WT_restedIL12_vs_individual_TP[["down_60"]] = union(down_all_WT_60minIL2_vs_WT_restedIL12,down_intron_WT_60minIL2_vs_WT_restedIL12)
WT_restedIL12_vs_individual_TP[["down_120"]] = union(down_all_WT_120minIL2_vs_WT_restedIL12,down_intron_WT_120minIL2_vs_WT_restedIL12)

KO_restedIL12_vs_individual_TP = list()
# up
up_all_KO_30minIL2_vs_KO_restedIL12 <- rownames(resALL[["norm_all"]][["KO_30minIL2_vs_KO_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_all_KO_60minIL2_vs_KO_restedIL12 <- rownames(resALL[["norm_all"]][["KO_60minIL2_vs_KO_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_all_KO_120minIL2_vs_KO_restedIL12 <- rownames(resALL[["norm_all"]][["KO_120minIL2_vs_KO_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
# down
down_all_KO_30minIL2_vs_KO_restedIL12 <- rownames(resALL[["norm_all"]][["KO_30minIL2_vs_KO_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_all_KO_60minIL2_vs_KO_restedIL12 <- rownames(resALL[["norm_all"]][["KO_60minIL2_vs_KO_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_all_KO_120minIL2_vs_KO_restedIL12 <- rownames(resALL[["norm_all"]][["KO_120minIL2_vs_KO_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))

# up
up_intron_KO_30minIL2_vs_KO_restedIL12 <- rownames(resALL[["norm_intron"]][["KO_30minIL2_vs_KO_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_KO_60minIL2_vs_KO_restedIL12 <- rownames(resALL[["norm_intron"]][["KO_60minIL2_vs_KO_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_KO_120minIL2_vs_KO_restedIL12 <- rownames(resALL[["norm_intron"]][["KO_120minIL2_vs_KO_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
# down
down_intron_KO_30minIL2_vs_KO_restedIL12 <- rownames(resALL[["norm_intron"]][["KO_30minIL2_vs_KO_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_KO_60minIL2_vs_KO_restedIL12 <- rownames(resALL[["norm_intron"]][["KO_60minIL2_vs_KO_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_KO_120minIL2_vs_KO_restedIL12 <- rownames(resALL[["norm_intron"]][["KO_120minIL2_vs_KO_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))

KO_restedIL12_vs_individual_TP[["up_30"]] = union(up_all_KO_30minIL2_vs_KO_restedIL12,up_intron_KO_30minIL2_vs_KO_restedIL12)
KO_restedIL12_vs_individual_TP[["up_60"]] = union(up_all_KO_60minIL2_vs_KO_restedIL12,up_intron_KO_60minIL2_vs_KO_restedIL12)
KO_restedIL12_vs_individual_TP[["up_120"]] = union(up_all_KO_120minIL2_vs_KO_restedIL12,up_intron_KO_120minIL2_vs_KO_restedIL12)

KO_restedIL12_vs_individual_TP[["down_30"]] = union(down_all_KO_30minIL2_vs_KO_restedIL12,down_intron_KO_30minIL2_vs_KO_restedIL12)
KO_restedIL12_vs_individual_TP[["down_60"]] = union(down_all_KO_60minIL2_vs_KO_restedIL12,down_intron_KO_60minIL2_vs_KO_restedIL12)
KO_restedIL12_vs_individual_TP[["down_120"]] = union(down_all_KO_120minIL2_vs_KO_restedIL12,down_intron_KO_120minIL2_vs_KO_restedIL12)


#up
s4 <- list( "WT_30up" = WT_restedIL12_vs_individual_TP[["up_30"]],
            "WT_60up" = WT_restedIL12_vs_individual_TP[["up_60"]],
            "WT_120up" = WT_restedIL12_vs_individual_TP[["up_120"]]
            )

pdf(paste0(out_dir_plots,"/09_WT_IL_up_changes_byTP.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

s4 <- list( "KO_30up" = KO_restedIL12_vs_individual_TP[["up_30"]],
            "KO_60up" = KO_restedIL12_vs_individual_TP[["up_60"]],
            "KO_120up" = KO_restedIL12_vs_individual_TP[["up_120"]]
            )

pdf(paste0(out_dir_plots,"/09_KO_IL_up_changes_byTP.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

#down
s4 <- list( "WT_30down" = WT_restedIL12_vs_individual_TP[["down_30"]],
            "WT_60down" = WT_restedIL12_vs_individual_TP[["down_60"]],
            "WT_120down" = WT_restedIL12_vs_individual_TP[["down_120"]]
            )

pdf(paste0(out_dir_plots,"/09_WT_IL_down_changes_byTP.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

s4 <- list( "KO_30down" = KO_restedIL12_vs_individual_TP[["down_30"]],
            "KO_60down" = KO_restedIL12_vs_individual_TP[["down_60"]],
            "KO_120down" = KO_restedIL12_vs_individual_TP[["down_120"]]
            )

pdf(paste0(out_dir_plots,"/09_KO_IL_down_changes_byTP.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

WT_up = unique(c(WT_restedIL12_vs_individual_TP[["up_30"]] ,WT_restedIL12_vs_individual_TP[["up_60"]] ,WT_restedIL12_vs_individual_TP[["up_120"]] ))
WT_down = unique(c(WT_restedIL12_vs_individual_TP[["down_30"]] ,WT_restedIL12_vs_individual_TP[["down_60"]] ,WT_restedIL12_vs_individual_TP[["down_120"]] ))

KO_up = unique(c(KO_restedIL12_vs_individual_TP[["up_30"]] ,KO_restedIL12_vs_individual_TP[["up_60"]] ,KO_restedIL12_vs_individual_TP[["up_120"]] ))
KO_down = unique(c(KO_restedIL12_vs_individual_TP[["down_30"]] ,KO_restedIL12_vs_individual_TP[["down_60"]] ,KO_restedIL12_vs_individual_TP[["down_120"]] ))

s4 <- list( "WT_up" =WT_up,
            "KO_up" = KO_up
            )

pdf(paste0(out_dir_plots,"/10_KO_vs_WT_anyUP.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

s4 <- list( "WT_down" = WT_down,
            "KO_down" = KO_down
            )

pdf(paste0(out_dir_plots,"/10_KO_vs_WT_anyDOWN.pdf"))
        p <- plot(euler(s4, shape = "ellipse"), quantities = TRUE)
        print(p)
dev.off()

# any up KO vs WT

KO_vs_WT_by_TP = list()

# up
up_all_KO_rested_vs_WT_rested <- rownames(resALL[["norm_all"]][["KO_restedIL12_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_all_KO_30minIL2_vs_WT_30minIL2 <- rownames(resALL[["norm_all"]][["KO_30minIL2_vs_WT_30minIL2"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_all_KO_60minIL2_vs_WT_60minIL2 <- rownames(resALL[["norm_all"]][["KO_60minIL2_vs_WT_60minIL2"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_all_KO_120minIL2_vs_WT_120minIL2 <- rownames(resALL[["norm_all"]][["KO_120minIL2_vs_WT_120minIL2"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))

# down
down_all_KO_rested_vs_WT_rested <- rownames(resALL[["norm_all"]][["KO_restedIL12_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_all_KO_30minIL2_vs_WT_30minIL2 <- rownames(resALL[["norm_all"]][["KO_30minIL2_vs_WT_30minIL2"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_all_KO_60minIL2_vs_WT_60minIL2 <- rownames(resALL[["norm_all"]][["KO_60minIL2_vs_WT_60minIL2"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_all_KO_120minIL2_vs_WT_120minIL2 <- rownames(resALL[["norm_all"]][["KO_120minIL2_vs_WT_120minIL2"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))

# up
up_intron_KO_rested_vs_WT_rested <- rownames(resALL[["norm_intron"]][["KO_restedIL12_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_KO_30minIL2_vs_WT_30minIL2 <- rownames(resALL[["norm_intron"]][["KO_30minIL2_vs_WT_30minIL2"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_KO_60minIL2_vs_WT_60minIL2 <- rownames(resALL[["norm_intron"]][["KO_60minIL2_vs_WT_60minIL2"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))
up_intron_KO_120minIL2_vs_WT_120minIL2 <- rownames(resALL[["norm_intron"]][["KO_120minIL2_vs_WT_120minIL2"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) > l2fc & as.numeric(as.character(padj)) <= pval))

# down
down_intron_KO_rested_vs_WT_rested <- rownames(resALL[["norm_intron"]][["KO_restedIL12_vs_WT_restedIL12"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_KO_30minIL2_vs_WT_30minIL2 <- rownames(resALL[["norm_intron"]][["KO_30minIL2_vs_WT_30minIL2"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_KO_60minIL2_vs_WT_60minIL2 <- rownames(resALL[["norm_intron"]][["KO_60minIL2_vs_WT_60minIL2"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))
down_intron_KO_120minIL2_vs_WT_120minIL2 <- rownames(resALL[["norm_intron"]][["KO_120minIL2_vs_WT_120minIL2"]] %>% dplyr::filter( as.numeric(as.character(log2FoldChange)) < -1*l2fc & as.numeric(as.character(padj)) <= pval))

KO_vs_WT_by_TP[["rested_bl_up"]] = union(up_all_KO_rested_vs_WT_rested,up_intron_KO_rested_vs_WT_rested)
KO_vs_WT_by_TP[["up_30"]] = union(up_all_KO_30minIL2_vs_WT_30minIL2,up_intron_KO_30minIL2_vs_WT_30minIL2)
KO_vs_WT_by_TP[["up_60"]] = union(up_all_KO_60minIL2_vs_WT_60minIL2,up_intron_KO_60minIL2_vs_WT_60minIL2)
KO_vs_WT_by_TP[["up_120"]] = union(up_all_KO_120minIL2_vs_WT_120minIL2,up_intron_KO_120minIL2_vs_WT_120minIL2)

KO_vs_WT_by_TP[["rested_bl_down"]] = union(down_all_KO_rested_vs_WT_rested,down_intron_KO_rested_vs_WT_rested)
KO_vs_WT_by_TP[["down_30"]] = union(down_all_KO_30minIL2_vs_WT_30minIL2,down_intron_KO_30minIL2_vs_WT_30minIL2)
KO_vs_WT_by_TP[["down_60"]] = union(down_all_KO_60minIL2_vs_WT_60minIL2,down_intron_KO_60minIL2_vs_WT_60minIL2)
KO_vs_WT_by_TP[["down_120"]] = union(down_all_KO_120minIL2_vs_WT_120minIL2,down_intron_KO_120minIL2_vs_WT_120minIL2)

## ANY:

ANY_change = unique(c( unlist(KO_vs_WT_by_TP),
                        unlist(KO_restedIL12_vs_individual_TP),
                        unlist(WT_restedIL12_vs_individual_TP) ))

mat_all = resALL[["norm_all"]][["Norm_expression"]]
mat_all = as.data.frame(mat_all[is.finite(rowSums(mat_all)),],stringsAsFactors=FALSE)
size = Ref_Anno_rn[rownames(mat_all),"width"] / 1000
mat_all = mat_all/size
mat_all[is.na(mat_all)] = 0
mat_all = mat_all[ANY_change ,c(   "WT_restedIL12","WT_30minIL2","WT_60minIL2","WT_120minIL2","KO_restedIL12","KO_30minIL2","KO_60minIL2","KO_120minIL2" )]
mat_all = as.data.frame(t(scale(t(mat_all),scale=TRUE,center=TRUE))[,colnames(mat_all)],stringsAsFactors=FALSE)

Plot = mat_all
sel_g = rep("ns",length( ANY_change ))
names(sel_g)= ANY_change

split_c = colnames(Plot)

split_c[grep("WT_restedIL12",split_c)] = "01_WT"
split_c[grepl("WT_30minIL2",split_c)] = "01_WT"
split_c[grep("WT_60minIL2",split_c)] = "01_WT"
split_c[grepl("WT_120minIL2",split_c)] = "01_WT"

split_c[grep("KO_restedIL12",split_c)] = "02_KO"
split_c[grepl("KO_30minIL2",split_c)] = "02_KO"
split_c[grep("KO_60minIL2",split_c)] = "02_KO"
split_c[grepl("KO_120minIL2",split_c)] = "02_KO"


sel_WT_anyIL2_change = sel_g

sel_WT_anyIL2_change[names(sel_WT_anyIL2_change) %in%     WT_down ] = "down"
sel_WT_anyIL2_change[names(sel_WT_anyIL2_change) %in%     WT_up ] = "up"

sel_KO_anyIL2_change = sel_g

sel_KO_anyIL2_change[names(sel_KO_anyIL2_change) %in%     KO_down ] = "down"
sel_KO_anyIL2_change[names(sel_KO_anyIL2_change) %in%     KO_up ] = "up"
    
sel_KO_vs_WT_diff_BaseLine = sel_g
sel_KO_vs_WT_diff_BaseLine[names(sel_KO_vs_WT_diff_BaseLine) %in%     KO_vs_WT_by_TP[["rested_bl_up"]] ] = "up"
sel_KO_vs_WT_diff_BaseLine[names(sel_KO_vs_WT_diff_BaseLine) %in%     KO_vs_WT_by_TP[["rested_bl_down"]]] = "down"

sel_KO_vs_WT_diff_IL2 = sel_g
sel_KO_vs_WT_diff_IL2[names(sel_KO_vs_WT_diff_IL2) %in%     unique(unlist(KO_vs_WT_by_TP[grep("up_",names(KO_vs_WT_by_TP))])) ] = "up"
sel_KO_vs_WT_diff_IL2[names(sel_KO_vs_WT_diff_IL2) %in%     unique(unlist(KO_vs_WT_by_TP[grep("down_",names(KO_vs_WT_by_TP))])) ] = "down"

ra= rowAnnotation(
                    sel_WT_anyIL2_change = sel_WT_anyIL2_change ,
                    sel_KO_anyIL2_change = sel_KO_anyIL2_change ,
                    sel_KO_vs_WT_diff_BaseLine = sel_KO_vs_WT_diff_BaseLine,
                    sel_KO_vs_WT_diff_IL2 = sel_KO_vs_WT_diff_IL2,

                    col = list( sel_WT_anyIL2_change = c("ns" = "white", "up" = "red", "down" = "blue"),
                                sel_KO_anyIL2_change = c("ns" = "white", "up" = "red", "down" = "blue"),
                                sel_KO_vs_WT_diff_BaseLine = c("ns" = "white", "up" = "red", "down" = "blue"),
                                sel_KO_vs_WT_diff_IL2 = c("ns" = "white", "up" = "red", "down" = "blue")  ),
                    Cyril_FCSinduced_sensitiveSRFlinked = anno_mark(at = which(rownames(Plot) %in% Cyril_FCSinduced_sensitiveSRFlinked), labels =rownames(Plot)[which(rownames(Plot) %in% Cyril_FCSinduced_sensitiveSRFlinked) ], labels_gp=gpar(fontsize = 1) ),
                    Gualdrini_SRFdirect = anno_mark(at = which(rownames(Plot) %in% Gualdrini_SRFdirect), labels =rownames(Plot)[which(rownames(Plot) %in% Gualdrini_SRFdirect) ], labels_gp=gpar(fontsize = 1) ),
                    annotation_name_gp = gpar(fontsize =3),
                    border=TRUE
                )



mat_clust = Plot
mat_clust[is.na(mat_clust)] = 0



# MFA - UMAP:
GG = c( length(colnames(mat_clust)[ grep("_restedIL12",colnames(mat_clust)) ]) ,
        length(colnames(mat_clust)[ grep("_30minIL2",colnames(mat_clust)) ]),
        length(colnames(mat_clust)[ grep("_60minIL2",colnames(mat_clust)) ]) ,
        length(colnames(mat_clust)[ grep("_120minIL2",colnames(mat_clust)) ]) )

MFA <- MFA(     cbind(  mat_clust[,grep("_restedIL12",colnames(mat_clust)) ],
                        mat_clust[,grep("_30minIL2",colnames(mat_clust)) ],
                        mat_clust[,grep("_60minIL2",colnames(mat_clust)) ],
                        mat_clust[,grep("_120minIL2",colnames(mat_clust)) ]) ,
                group = GG,type = rep("c",length(GG)),
                ncp=ncol(mat_clust),
                name.group=c("_restedIL12","_30minIL2","_60minIL2","_120minIL2"),
                graph=FALSE)

indiv = MFA$ind$coord

rpca = FactoMineR::PCA(mat_clust, graph = FALSE,scale.unit = FALSE,ncp = ncol(mat_clust))
indiv = rpca$ind$coord[,1:(dim(rpca$eig[rpca$eig[,3]<=99,]))[1] ]

k = round(sqrt( dim(indiv)[1] * dim(indiv)[2]  ) )

knn.real = get.knn(as.matrix(indiv), k = k, algorithm="brute")
knn.real = data.frame(from = rep(1:nrow(knn.real$nn.index), k), to = as.vector(knn.real$nn.index), weight = 1/(1 + as.vector(knn.real$nn.dist)))
nw.real = graph_from_data_frame(knn.real, directed = FALSE)
nw.real = simplify(nw.real)
lc.real = cluster_louvain(nw.real)
louvain = as.factor(membership(lc.real)) 


cl = paste0("Cluster_",louvain)
names(cl) = rownames(indiv)
cln=cl

HT = Heatmap(     Plot,
                        width = unit(20, "mm"),
                        row_split = paste0("Cluster_",cl) , 
                        column_split = split_c,
                        column_gap = unit(1, "mm"),
                        row_gap = unit(1, "mm"),
                        cluster_column_slices = FALSE,
                        cluster_columns=FALSE, 
                        right_annotation = ra,
                        show_row_names = FALSE,
                        show_column_names = TRUE,
                        column_title_gp = gpar(fontsize = 3),
                        row_title_gp = gpar(fontsize = 3),
                        heatmap_legend_param = list(title = "z-scored ATACsignal",direction = "horizontal"),
                        column_names_gp = gpar(fontsize = 2),
                        row_dend_width = unit(2, "mm"),
                        show_row_dend = TRUE,
                        use_raster = TRUE,
                        row_dend_reorder = TRUE,
                        column_dend_reorder = FALSE,
                        na_col = "white",
                        border=TRUE,
                        col = colorRamp2( seq(-2, 2, by = 0.1), cividis(length(seq(-2, 2, by = 0.1))))
                        )
pdf(paste0(out_dir_plots,"/11_Heat_all_intron.pdf"),useDingbats = FALSE)
    draw(HT, heatmap_legend_side="right")
dev.off()

Plot$Clusters = NA
Plot[names(cln) , "Clusters" ] = cln
write.table(Plot,paste0(out_dir_plots,"/11_Heat_Clusters.txt"),sep="\t",col.names=NA)
Plot_11 = Plot


cl = cln[names(cln) %in% Universe]
list_comp = split(names(cl),cl)
check = unique(unlist(list_comp))
Universe = union(Universe,check)

# for each list assess the overlap:
# HG + JI:

msig_plots = paste0(out_dir_plots,"/msig_HJ_IL2_TC/")
unlink(msig_plots, recursive=TRUE)
dir.create(msig_plots)

# for each list assess the overlap:
# HG + JI:
HGJc = lapply( 1:length(MSIG_list_symbol_entrez_mouse), function(j)  {

     l1 = lapply(1:length(list_comp),function(lv) {
        l2 = lapply( 1:length(MSIG_list_symbol_entrez_mouse[[j]]), function(y)  {
            set1=list_comp[[lv]] #gsub(".*_","",list_comp[[lv]])
            set2=MSIG_list_symbol_entrez_mouse[[j]][[y]]#gsub(".*_","",MSIG_list_symbol_entrez_mouse[[j]][[y]])

            HGu=HG_fun_upper(set1,set2,Universe)
            HGl=HG_fun_lower(set1,set2,Universe)
            HG = min(c(HGu,HGl))

            JI= JI_fun(set1,set2)

            Genes = paste0(intersect(set1,set2),collapse=";")

            j_nm<-c()
            for (nm_rep in 1:500){
                A_<-sample(Universe,length(set1))
                p_rand = JI_fun(A_,set2)
                j_nm<-c(j_nm,p_rand)
            }

            SCALED_JI<-(JI-mean(j_nm))/sd(j_nm)
            if(is.nan(SCALED_JI)){SCALED_JI=0}
            return(c( HG=HG,HGu=HGu,HGl=HGl,JI= JI, SCALED_JI = SCALED_JI,Genes=Genes ))

        })
        names(l2) = names(MSIG_list_symbol_entrez_mouse[[j]])
        l2 = as.data.frame(do.call(rbind,l2),stringsAsFactors=FALSE)
        HG = as.numeric(as.character(l2$HG)) #
        names(HG) = rownames(l2)        
        HGu = as.numeric(as.character(l2$HGu)) #
        names(HGu) = rownames(l2)
        HGl = as.numeric(as.character(l2$HGl)) #
        names(HGl) = rownames(l2)
        JI = as.numeric(as.character(l2$JI))
        names(JI) = rownames(l2)
        SCALED_JI = as.numeric(as.character(l2$SCALED_JI))
        names(SCALED_JI) = rownames(l2)
        Genes = l2$Genes
        names(Genes) = rownames(l2)
        return(list(HG =HG, HGu = HGu,HGl=HGl,JI=JI,SCALED_JI=SCALED_JI,Genes=Genes))
    })

    names(l1) = names(list_comp)
    HG = as.data.frame(do.call(cbind,lapply(l1, `[[`, 1)),stringsAsFactors=FALSE)    
    HGu = as.data.frame(do.call(cbind,lapply(l1, `[[`, 2)),stringsAsFactors=FALSE)
    HGl = as.data.frame(do.call(cbind,lapply(l1, `[[`, 3)),stringsAsFactors=FALSE)
    JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 4)),stringsAsFactors=FALSE)
    SCALED_JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 5)),stringsAsFactors=FALSE)
    Genes = as.data.frame(do.call(cbind,lapply(l1, `[[`, 6)),stringsAsFactors=FALSE)

    res_l = list(   
                    HG = as.data.frame(do.call(cbind,lapply(l1, `[[`, 1)),stringsAsFactors=FALSE)    ,
                    HGu = as.data.frame(do.call(cbind,lapply(l1, `[[`, 2)),stringsAsFactors=FALSE),
                    HGl = as.data.frame(do.call(cbind,lapply(l1, `[[`, 3)),stringsAsFactors=FALSE),
                    JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 4)),stringsAsFactors=FALSE),
                    SCALED_JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 5)),stringsAsFactors=FALSE),
                    FDRu = apply(HGu,2,function(x){p.adjust(x,method='BH')}),
                    FDRl = apply(HGl,2,function(x){p.adjust(x,method='BH')}),
                    Genes = as.data.frame(do.call(cbind,lapply(l1, `[[`, 6)),stringsAsFactors=FALSE)
                    )

    output_file <- paste0(msig_plots,"/",names(MSIG_list_symbol_entrez_mouse)[j], "_MSIG_Enrichment_Results.xlsx")
    write.xlsx( res_l, output_file ,row.names = TRUE )

    return(list(HG=HG,HGu = HGu,HGl=HGl,JI=JI,SCALED_JI=SCALED_JI,Genes=Genes))

})

names(HGJc) = names(MSIG_list_symbol_entrez_mouse)


for(xx in names(MSIG_list_symbol_entrez_mouse)){

    HGT = HGJc[[xx]][["HG"]]
    HGT = apply(HGT,2,function(x){p.adjust(x,method='BH')})
    HGT = as.data.frame(HGT,stringsAsFactors=FALSE)
    HGT=-log10(HGT)
    max_replace = max(c(max(unlist(HGT)[is.finite(unlist(HGT))])))
    HGT[] <- lapply(HGT, function(i) if(is.numeric(i)) ifelse(is.infinite(i), max_replace, i) else i)
    if(length(ww)!=0){
    ww = which(apply(HGT, MARGIN=c(1), max) >= 2   )
    HGT = HGT[ww,,drop=FALSE]
    compress = 3
    HGT[HGT>=compress] = compress
    HGTrad = sqrt(HGT/pi)
    HGTrad = HGTrad/max(HGTrad)

    SCALED_JI = HGJc[[xx]][["SCALED_JI"]]

    rr = intersect(rownames(SCALED_JI),rownames(HGTrad))
    cc = intersect(colnames(SCALED_JI),colnames(HGTrad))

    SCALED_JI = SCALED_JI[rr,cc,drop=FALSE]
    HGTrad = HGTrad[rr,cc,drop=FALSE]

    if(dim(HGTrad)[1] !=0 & dim(HGTrad)[2] !=0){

        col_fun = colorRamp2(c(-1*max(abs(SCALED_JI)),-1*max(abs(SCALED_JI))/2,0,max(abs(SCALED_JI))/2,max(abs(SCALED_JI))),c("#034961","#51D06D", "white","#FF9E39","#C51300")  )
        
        htord = Heatmap(    SCALED_JI, 
                                cluster_rows = TRUE, 
                                row_dend_reorder = TRUE,
                                cluster_columns = TRUE,
                                column_dend_reorder=TRUE
                            )
        cO = colnames(SCALED_JI)[unlist(column_order(htord))]
        rO = rownames(SCALED_JI)[unlist(row_order(htord))]
        SCALED_JI = SCALED_JI[rO,cO,drop=FALSE]
        HGTrad = HGTrad[rownames(SCALED_JI),colnames(SCALED_JI),drop=FALSE]
        dds = 20/max(nrow(SCALED_JI),ncol(SCALED_JI))
               HR1 = Heatmap(          SCALED_JI,

                                column_title = paste0("Fraction: radius min: ",round(min(HGJc[[xx]][["HG"]]),4)," radius max: ",round(max(HGJc[[xx]][["HG"]]),4)),
                                width = ncol(SCALED_JI)*unit(dds, "cm"), 
                                height = nrow(SCALED_JI)*unit(dds, "cm"),
                                #rect_gp = gpar(type = "none"), 
                                #cell_fun = function(j, i, x, y, width, height, fill) {
                                #    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = NA, fill = NA))
                                #    grid.circle(x = x, y = y, r = abs(HGTrad[i, j])/2.5 * unit(dds, "cm"), 
                                #    gp = gpar(fill = col_fun(SCALED_JI[i, j]), col = NA))
                                #}, 
                                show_row_names = TRUE,
                                show_column_names = TRUE,
                                column_title_gp = gpar(fontsize = 2),
                                row_title_gp = gpar(fontsize = 1),
                                row_names_gp = gpar(fontsize = 1),
                                heatmap_legend_param = list(title = "Calls",direction = "horizontal"),
                                column_names_gp = gpar(fontsize = 2),
                                row_dend_width = unit(2, "mm"),
                                border=TRUE,
                                col = col_fun ,
                                cluster_rows = TRUE, 
                                cluster_columns = TRUE,
                                row_dend_reorder = TRUE,
                                column_dend_reorder=TRUE

                                )
        pdf(paste0(msig_plots,"/", gsub(" |:","_",xx) ,"_enrichment.pdf"),useDingbats = FALSE , width=10   , height=10)
                draw(HR1)
        dev.off()  
    }
}
}

#####

ref_list = split(rownames(Plot_11),Plot_11$Clusters)
Plot_08_list = split(rownames(Plot_08),Plot_08$Clusters)
Plot_04_list = split(rownames(Plot_04),Plot_04$Clusters)

l1 = lapply(1:length(Plot_04_list),function(lv) {
    l2 = lapply( 1:length(ref_list), function(y)  {
        set1=Plot_04_list[[lv]] #gsub(".*_","",list_comp[[lv]])
        set2=ref_list[[y]]#gsub(".*_","",MSIG_list_symbol_entrez_mouse[[j]][[y]])

        HGu=HG_fun_upper(set1,set2,Universe)
        HGl=HG_fun_lower(set1,set2,Universe)
        HG = min(c(HGu,HGl))

        JI= JI_fun(set1,set2)

        j_nm<-c()
        for (nm_rep in 1:500){
            A_<-sample(Universe,length(set1))
            p_rand = JI_fun(A_,set2)
            j_nm<-c(j_nm,p_rand)
        }

        SCALED_JI<-(JI-mean(j_nm))/sd(j_nm)
        if(is.nan(SCALED_JI)){SCALED_JI=0}
        return(c( HG=HG,HGu=HGu,HGl=HGl,JI= JI, SCALED_JI = SCALED_JI ))

    })
    names(l2) = names(ref_list)
    l2 = as.data.frame(do.call(rbind,l2),stringsAsFactors=FALSE)
    HG =l2$HG #
    names(HG) = rownames(l2)        
    HGu =l2$HGu #
    names(HGu) = rownames(l2)
    HGl =l2$HGl #
    names(HGl) = rownames(l2)
    JI = l2$JI
    names(JI) = rownames(l2)
    SCALED_JI = l2$SCALED_JI
    names(SCALED_JI) = rownames(l2)
    return(list(HG =HG, HGu = HGu,HGl=HGl,JI=JI,SCALED_JI=SCALED_JI))
})

names(l1) = names(Plot_04_list)
HG = as.data.frame(do.call(cbind,lapply(l1, `[[`, 1)),stringsAsFactors=FALSE)    
HGu = as.data.frame(do.call(cbind,lapply(l1, `[[`, 2)),stringsAsFactors=FALSE)
HGl = as.data.frame(do.call(cbind,lapply(l1, `[[`, 3)),stringsAsFactors=FALSE)
JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 4)),stringsAsFactors=FALSE)
SCALED_JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 5)),stringsAsFactors=FALSE)

res_l = list(   
                HG = as.data.frame(do.call(cbind,lapply(l1, `[[`, 1)),stringsAsFactors=FALSE)    ,
                HGu = as.data.frame(do.call(cbind,lapply(l1, `[[`, 2)),stringsAsFactors=FALSE),
                HGl = as.data.frame(do.call(cbind,lapply(l1, `[[`, 3)),stringsAsFactors=FALSE),
                JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 4)),stringsAsFactors=FALSE),
                SCALED_JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 5)),stringsAsFactors=FALSE),
                FDRu = apply(HGu,2,function(x){p.adjust(x,method='BH')}),
                FDRl = apply(HGl,2,function(x){p.adjust(x,method='BH')})
                )

output_file <- paste0(paste0(out_dir_plots,"/msig_HJ_NAIVE_TCR24h/"),"/Enrichment_IL2Clusters.xlsx")
write.xlsx( res_l, output_file ,row.names = TRUE )


SCALED_JI = res_l[["SCALED_JI"]]
col_fun = colorRamp2(c( -1*max(abs(SCALED_JI))/2,
                        -1*max(abs(SCALED_JI))/4,
                        0,
                        max(abs(SCALED_JI))/4,
                        max(abs(SCALED_JI))/2),
                        c("#034961","#51D06D", "white","#FF9E39","#C51300")  )
htord = Heatmap(    SCALED_JI, 
                        cluster_rows = TRUE, 
                        row_dend_reorder = TRUE,
                        cluster_columns = TRUE,
                        column_dend_reorder=TRUE
                    )
cO = colnames(SCALED_JI)[unlist(column_order(htord))]
rO = rownames(SCALED_JI)[unlist(row_order(htord))]
SCALED_JI = SCALED_JI[rO,cO,drop=FALSE]

dds = 20/max(nrow(SCALED_JI),ncol(SCALED_JI))
HR1 = Heatmap(          SCALED_JI,
            width = ncol(SCALED_JI)*unit(dds, "cm"), 
            height = nrow(SCALED_JI)*unit(dds, "cm"),
            show_row_names = TRUE,
            show_column_names = TRUE,
            column_title_gp = gpar(fontsize = 2),
            row_title_gp = gpar(fontsize = 1),
            row_names_gp = gpar(fontsize = 1),
            heatmap_legend_param = list(title = "Calls",direction = "horizontal"),
            column_names_gp = gpar(fontsize = 2),
            row_dend_width = unit(2, "mm"),
            border=TRUE,
            col = col_fun ,
            cluster_rows = TRUE, 
            cluster_columns = TRUE,
            row_dend_reorder = TRUE,
            column_dend_reorder=TRUE

            )
pdf(paste0(out_dir_plots,"/Clusters_04_to_11_enrichment.pdf"),useDingbats = FALSE , width=10   , height=10)
        draw(HR1)
dev.off() 


l1 = lapply(1:length(Plot_08_list),function(lv) {
    l2 = lapply( 1:length(ref_list), function(y)  {
        set1=Plot_08_list[[lv]] #gsub(".*_","",list_comp[[lv]])
        set2=ref_list[[y]]#gsub(".*_","",MSIG_list_symbol_entrez_mouse[[j]][[y]])

        HGu=HG_fun_upper(set1,set2,Universe)
        HGl=HG_fun_lower(set1,set2,Universe)
        HG = min(c(HGu,HGl))

        JI= JI_fun(set1,set2)

        j_nm<-c()
        for (nm_rep in 1:500){
            A_<-sample(Universe,length(set1))
            p_rand = JI_fun(A_,set2)
            j_nm<-c(j_nm,p_rand)
        }

        SCALED_JI<-(JI-mean(j_nm))/sd(j_nm)
        if(is.nan(SCALED_JI)){SCALED_JI=0}
        return(c( HG=HG,HGu=HGu,HGl=HGl,JI= JI, SCALED_JI = SCALED_JI ))

    })
    names(l2) = names(ref_list)
    l2 = as.data.frame(do.call(rbind,l2),stringsAsFactors=FALSE)
    HG =l2$HG #
    names(HG) = rownames(l2)        
    HGu =l2$HGu #
    names(HGu) = rownames(l2)
    HGl =l2$HGl #
    names(HGl) = rownames(l2)
    JI = l2$JI
    names(JI) = rownames(l2)
    SCALED_JI = l2$SCALED_JI
    names(SCALED_JI) = rownames(l2)
    return(list(HG =HG, HGu = HGu,HGl=HGl,JI=JI,SCALED_JI=SCALED_JI))
})

names(l1) = names(Plot_08_list)
HG = as.data.frame(do.call(cbind,lapply(l1, `[[`, 1)),stringsAsFactors=FALSE)    
HGu = as.data.frame(do.call(cbind,lapply(l1, `[[`, 2)),stringsAsFactors=FALSE)
HGl = as.data.frame(do.call(cbind,lapply(l1, `[[`, 3)),stringsAsFactors=FALSE)
JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 4)),stringsAsFactors=FALSE)
SCALED_JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 5)),stringsAsFactors=FALSE)

res_l = list(   
                HG = as.data.frame(do.call(cbind,lapply(l1, `[[`, 1)),stringsAsFactors=FALSE)    ,
                HGu = as.data.frame(do.call(cbind,lapply(l1, `[[`, 2)),stringsAsFactors=FALSE),
                HGl = as.data.frame(do.call(cbind,lapply(l1, `[[`, 3)),stringsAsFactors=FALSE),
                JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 4)),stringsAsFactors=FALSE),
                SCALED_JI = as.data.frame(do.call(cbind,lapply(l1, `[[`, 5)),stringsAsFactors=FALSE),
                FDRu = apply(HGu,2,function(x){p.adjust(x,method='BH')}),
                FDRl = apply(HGl,2,function(x){p.adjust(x,method='BH')})
                )

output_file <- paste0(paste0(out_dir_plots,"/msig_HJ_TCR24h_RESTED/"),"/Enrichment_IL2Clusters.xlsx")
write.xlsx( res_l, output_file ,row.names = TRUE )



SCALED_JI = res_l[["SCALED_JI"]]
col_fun = colorRamp2(c( -1*max(abs(SCALED_JI))/2,
                        -1*max(abs(SCALED_JI))/4,
                        0,
                        max(abs(SCALED_JI))/4,
                        max(abs(SCALED_JI))/2),
                        c("#034961","#51D06D", "white","#FF9E39","#C51300")  )
htord = Heatmap(    SCALED_JI, 
                        cluster_rows = TRUE, 
                        row_dend_reorder = TRUE,
                        cluster_columns = TRUE,
                        column_dend_reorder=TRUE
                    )
cO = colnames(SCALED_JI)[unlist(column_order(htord))]
rO = rownames(SCALED_JI)[unlist(row_order(htord))]
SCALED_JI = SCALED_JI[rO,cO,drop=FALSE]

dds = 20/max(nrow(SCALED_JI),ncol(SCALED_JI))
HR1 = Heatmap(          SCALED_JI,
            width = ncol(SCALED_JI)*unit(dds, "cm"), 
            height = nrow(SCALED_JI)*unit(dds, "cm"),
            show_row_names = TRUE,
            show_column_names = TRUE,
            column_title_gp = gpar(fontsize = 2),
            row_title_gp = gpar(fontsize = 1),
            row_names_gp = gpar(fontsize = 1),
            heatmap_legend_param = list(title = "Calls",direction = "horizontal"),
            column_names_gp = gpar(fontsize = 2),
            row_dend_width = unit(2, "mm"),
            border=TRUE,
            col = col_fun ,
            cluster_rows = TRUE, 
            cluster_columns = TRUE,
            row_dend_reorder = TRUE,
            column_dend_reorder=TRUE

            )
pdf(paste0(out_dir_plots,"/Clusters_08_to_11_enrichment.pdf"),useDingbats = FALSE , width=10   , height=10)
        draw(HR1)
dev.off() 