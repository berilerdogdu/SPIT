options(repos = c(CRAN = "https://cloud.r-project.org"))
setwd(directory_path)

import_infreps <- function(pheno, tx2gene, quant_type, quant_path, filtered_tx_counts){
    targets <- read_tsv(pheno)
    path <- file.path(quant_path, targets$id, "quant.sf")
    if(all(file.exists(path))==FALSE){
        print("Quantification files not found, please make sure to specify the parent directory correctly.", flush=TRUE)
        return(NULL)
    }
    names(path) <- targets$id
    tx2gene_df <- read_tsv(tx2gene)
    Txi_tx <- tximport(path,
    type = quant_type,
    tx2gene = tx2gene_df,
    txOut = TRUE,
    countsFromAbundance = "dtuScaledTPM",
    ignoreTxVersion = FALSE)
    
    filtered_txs <- rownames(read.table(filtered_tx_counts, header = TRUE, sep = "\t", row.names = 1))
    pb <- progress_bar$new(total = length(names(Txi_tx$infReps)), format = "[:bar] :percent :elapsed")
    for (mat_name in names(Txi_tx$infReps)) {
        mat <- Txi_tx$infReps[[mat_name]]
        rownames(mat) <- rownames(Txi_tx$counts)
        selected_rows <- rownames(mat) %in% filtered_txs
        mat <- mat[selected_rows, ]
        file_name <- paste("tximport_infReps_sample_", toString(mat_name), ".txt", sep = "")
        write.table(mat, file.path("infReps", file_name), sep = '\t')
        pb$tick()
    }
    pb$terminate()
}
