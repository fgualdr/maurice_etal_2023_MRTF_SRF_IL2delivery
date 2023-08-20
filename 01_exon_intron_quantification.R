
library( renv )
library( AnnotationDbi )
library( GenomicFeatures )
library( GenomicAlignments )
library( Rsamtools )
library( BiocParallel )
library( openxlsx )
library( tidyverse )
library( yaml )

Ref_Txdb <- loadDb( "MM10_Txdb_mergedTx.sqlite" )
babs_file <- file.path( "..", "..", "..", ".babs" )
babs <- read_yaml( babs_file )[[ 1 ]] %>%
    map( paste, collapse = ",", sep = "" )
babs_df <- data.frame( label = names( babs ),
                       value = unlist( babs ) )
# The GFF/sqlite is in UCSC format

TxbyGene <- transcriptsBy(Ref_Txdb,
                          by="gene")

ExonByGene <- exonsBy(Ref_Txdb,
                      by="gene")

IntronByGene <- intronsByTranscript(Ref_Txdb,use.names=TRUE)

all_bams <- list.files( file.path( "..", "..", "nfcore", "mm10_ucsc_igenomes", "star_rsem" ),
                        pattern = ".*.markdup.sorted.bam$",
                        full.names = TRUE )

bamfiles <- BamFileList(all_bams)

register(MulticoreParam(workers=16))

se_all <- summarizeOverlaps(features=TxbyGene,
    reads=bamfiles,
    mode="Union",
    singleEnd=TRUE,
    ignore.strand=TRUE, 
    preprocess.reads=invertStrand,
    ## The inversion depends on the type of library made for the one recieved from you last time I needed to invert as it was stranded but complementary.
    inter.feature=FALSE)

dat_all <- as.data.frame(assay(se_all)) %>%
    rownames_to_column( var = "gene_id" )

se_int <- summarizeOverlaps(features=IntronByGene,
    reads=bamfiles,
    mode="IntersectionNotEmpty",
    singleEnd=TRUE,
    ignore.strand=TRUE,
    preprocess.reads=invertStrand,
    inter.feature=FALSE)

dat_int <- as.data.frame(assay(se_int)) %>%
    rownames_to_column( var = "gene_id" )

se_exonic <- summarizeOverlaps(features=ExonByGene,
    reads=bamfiles,
    mode="IntersectionStrict",
    singleEnd=TRUE,
    ignore.strand=TRUE,
    preprocess.reads=invertStrand,
    inter.feature=FALSE)

dat_ex <- as.data.frame(assay(se_exonic)) %>%
    rownames_to_column( var = "gene_id" )

seqlevels_tx <- data.frame( seqlevels = seqlevels( TxbyGene ) )

seqlevels_genome <- read.table( file = file.path( "..", "..", "docs", "mm10_ucsc_igenomes.fa.fai" ) ) %>%
    rename( seqlevels = V1 ) %>%
    dplyr::select( seqlevels )

seqlevels_comp <- list( seqlevels_tx = seqlevels_tx,
                       seqlevels_genome = seqlevels_genome ) %>%
    bind_rows( ) %>%
    distinct( ) %>%
    mutate( seqlevels_tx = seqlevels %in% seqlevels_tx$seqlevels ) %>%
    mutate( seqlevels_genome = seqlevels %in% seqlevels_genome$seqlevels )

output <- list( all = dat_all,
               intron = dat_int,
               exon = dat_ex,
               seqlevels = seqlevels_comp,
               babs_info = babs_df )
               
output_file <- paste0( "quantifcation_exon_intron_mm10_ucsc_igenome_",
                      sub( "@crick.ac.uk", "", babs$Scientist ),
                      "_", babs$Lims, ".xlsx" )
write.xlsx( output, output_file )
 
