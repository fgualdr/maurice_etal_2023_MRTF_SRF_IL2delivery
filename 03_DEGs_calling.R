library(openxlsx)
library(DESeq2)
library(dplyr)
library(ImpulseDE2)
library(zoo)

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


anno_path = "MM10_merged_annotation.txt"
Ref_Anno <- read.delim( anno_path,sep="\t" ,row.names=1)
Ref_Anno$id = paste0(rownames(Ref_Anno),"_",Ref_Anno$Entrez.Gene.ID)

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
resALL = list()
for(nn in names(list_mat)){
    cat(nn,"\n")
    table_out_dir = paste0(out_dir,"/",nn,"_DESEQ2/")
    dir.create(table_out_dir)
    mat = list_mat[[nn]]
    rownames(mat) = Ref_Anno[rownames(mat),"id"]
    mat_norm = sweep(mat[,Norm_mat$Sample_ID], 2, Norm_mat[,"scaling"], "*")
    Ref_Anno_rn = Ref_Anno
    rownames(Ref_Anno_rn) = Ref_Anno_rn$id
    if(nn == "norm_all"){
        size = Ref_Anno_rn[rownames(mat_norm),"width"] / 1000
    }else{
        if(nn == "norm_exon" ){
            size = Ref_Anno_rn[rownames(mat_norm),"cDNA_LENGTH"] / 1000
        }else{
            if(nn == "norm_intron"){
                size = (Ref_Anno_rn[rownames(mat_norm),"width"] / 1000)-(Ref_Anno_rn[rownames(mat_norm),"cDNA_LENGTH"] / 1000)
            }else{stopp("ERROR!")}
        }
    }
    limit_id = (mat_norm[,rownames(deg_design)]/size)
    AV = as.data.frame(matrix(NA,ncol=length(unique(deg_design$Sample_Condition)),nrow= nrow(limit_id) ),stringsAsFactors=FALSE)
    colnames(AV) = unique(deg_design$Sample_Condition)
    rownames(AV) = rownames(limit_id)
    for(cc in unique(deg_design$Sample_Condition)){
        sel = rownames(deg_design)[deg_design$Sample_Condition %in% cc]
        w1 = which(colnames(limit_id) %in% sel)
        AV[,cc] = rowMeans(limit_id[,w1])
    }
    signal_rna = as.numeric(as.character(log(as.matrix(AV))))
    signal_rna = signal_rna[!is.na(signal_rna)]
    signal_rna = signal_rna[is.finite(signal_rna)]
    skew_mod <-  sn::selm(signal_rna ~ 1)
    modeS = sn::modeSECdistr(skew_mod@param$dp,family="SN")
    limit_id[limit_id==0] = min(limit_id[limit_id!=0])
    TxAll = (log( limit_id )- modeS ) /skew_mod@param$cp["s.d."]
    rn = lapply(split_tibble(deg_design,"Sample_Condition"),function(x){
        return(rows_in_lim(TxAll[,rownames(x)],0.01,5))
    })
    rn =unique(unlist(rn))
    mat_sel = round(mat_norm[rn,],0)
    w = apply(mat_sel,1,mean)
    w = which(w > 100 & w < 1.5e+05)
    mat_sel = mat_sel[names(w),]

    # Comparisons:
    
    formula = ~ Sample_Replicate + Sample_Condition
    cond1_revel = "WT_naive"

    dds_ex <- DESeqDataSetFromMatrix(countData = mat_sel[,rownames(deg_design)], colData = deg_design, design = formula)
    dds_ex <- estimateSizeFactors(dds_ex)
    colData(dds_ex)$sizeFactor <- 1
    dds_ex$Sample_Condition = relevel(   dds_ex$Sample_Condition, cond1_revel)

    dds_ex <- estimateDispersions(dds_ex,fitType="local",maxit=300000) # parametric local mean
    dds_ex <- nbinomWaldTest(dds_ex,maxit = 300000)       
    pdf(paste0(table_out_dir,"/MAplot.pdf"))
        plotDispEsts(dds_ex)
    dev.off()
    RR = resultsNames(dds_ex)
    RRdat = as.data.frame(cbind(RR,rep(0,length(RR))))

    cc = as.character( unique(deg_design$Sample_Condition) )

    ll = list(  c1=c("WT_TCR24h","WT_naive"),
                c1=c("KO_TCR24h","KO_naive"),
                c1=c("WT_restedIL12","WT_TCR24h"),
                c1=c("KO_restedIL12","KO_TCR24h"),

                c1=c("WT_30minIL2","WT_restedIL12"),
                c1=c("WT_60minIL2","WT_restedIL12"),
                c1=c("WT_120minIL2","WT_restedIL12"),

                c1=c("KO_30minIL2","KO_restedIL12"),
                c1=c("KO_60minIL2","KO_restedIL12"),
                c1=c("KO_120minIL2","KO_restedIL12"),

                c1=c("KO_naive","WT_naive"),
                c1=c("KO_TCR24h","WT_TCR24h"),
                c1=c("KO_restedIL12","WT_restedIL12"),
                c1=c("KO_30minIL2","WT_30minIL2"),
                c1=c("KO_60minIL2","WT_60minIL2"),
                c1=c("KO_120minIL2","WT_120minIL2")

            )


    contrasts = do.call(rbind,ll)

    res_l = list( Norm_expression=AV )
    cat("Extract by contrast\n")
    for(rn in 1:nrow(contrasts)){

        res = results(dds_ex, contrast=c("Sample_Condition",contrasts[rn,1],contrasts[rn,2]))
        nm <- paste0(table_out_dir,"/",contrasts[rn,1],"_vs_",contrasts[rn,2],".txt")
        pdf(paste0(table_out_dir,"/",contrasts[rn,1],"_vs_",contrasts[rn,2],".pdf"))
                plotMA(res,alpha = 0.01, main = "", xlab = "mean of normalized counts", MLE = FALSE)
        dev.off()
        res = as.data.frame(res,stringsAsFactors=FALSE)
        cn_samp = deg_design[grep(paste0(contrasts[rn,2],"|",contrasts[rn,1]),deg_design$Sample_Condition), ]
        AV_sel = AV[,as.character(unique(cn_samp$Sample_Condition))]
        res = cbind(AV_sel[rownames(AV_sel),],res[rownames(AV_sel),])
        write.table(res,file=nm,sep="\t",col.names=NA)

        res_l[[ paste0(contrasts[rn,1],"_vs_",contrasts[rn,2] )  ]] = res

    }

    output_file <- paste0(table_out_dir,"/Master_Results.xlsx")
    write.xlsx( res_l, output_file ,row.names = TRUE )

    resALL[[nn]] = res_l

}

# Save RDat

save(resALL,file=paste0(out_dir,"/ListDat.RData"))


# Plots:
