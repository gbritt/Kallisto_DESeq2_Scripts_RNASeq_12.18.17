load_transcripts2gene_mappings = function(){
    library('biomaRt')
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                             dataset = "scerevisiae_gene_ensembl",
                             host = 'ensembl.org')
    Sys.sleep(2)
    t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                         "external_gene_name", "description"), mart = mart)
    t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                         ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
}
