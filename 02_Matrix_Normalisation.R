library( BiocParallel )
library( openxlsx )


anno_path = "MM10_merged_annotation.txt"
Ref_Anno <- read.delim( anno_path,sep="\t" )

### Diane_excel : 
out_dir = "/maurice_etal_2023_MRTF_SRF_IL2delivery/"
deg_design = read.delim("/maurice_etal_2023_MRTF_SRF_IL2delivery/DEG_DESIGN",sep="\t")
Diane_excel = "/maurice_etal_2023_MRTF_SRF_IL2delivery/quantifcation_exon_intron_mm10_ucsc_igenomes_diane.maurice_asf-RN22062.xlsx"
# Make individual data_tables all, intron, exons

all = read.xlsx(Diane_excel, "all",rowNames=TRUE,colNames=TRUE)
exon = read.xlsx(Diane_excel, "exon",rowNames=TRUE,colNames=TRUE)
intron = all - exon[rownames(all),colnames(all)]

## 
library(GeneralNormalizer)
param <- BiocParallel::MulticoreParam(workers=2,progressbar = TRUE)
Result = RunNorm(all,deg_design,fix_reference="random",row_name_index=1,saving_path=out_dir,n_pop=1,BiocParam=param)

# QC_plot
QC_plot(Result,deg_design,saving_path=out_dir)

# Results is an gnmsn_obj that can be piped into the QC 
# Make normalised all, exon, intron
scaling_factors = slot(Result, "scaling_factors")
norm_all = sweep(all[,scaling_factors$Sample_ID], 2, scaling_factors[,"scaling"], "*")
norm_exon = sweep(exon[,scaling_factors$Sample_ID], 2, scaling_factors[,"scaling"], "*")
norm_intron = sweep(intron[,scaling_factors$Sample_ID], 2, scaling_factors[,"scaling"], "*")
output <- list( all_norm = norm_all,
               intron_norm = norm_intron,
               exon_norm = norm_exon,
               scaling_factors = scaling_factors
                )
output_file <- paste0(out_dir, "Normalised_quantifcation_all_exon_intron_mm10_fgualdrini.xlsx")
write.xlsx( output, output_file )

# 