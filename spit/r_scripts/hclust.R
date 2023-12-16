options(repos = c(CRAN = "https://cloud.r-project.org"))
setwd(directory_path)

perform_hclust <- function(pheno, include_shared_dtu, col_pallete="BuGn", cov){
    controlled_file_path <- "controlled_spit_cluster_matrix.txt"
    non_controlled_file_path <- "spit_cluster_matrix.txt"
    if (file.exists(controlled_file_path)) {
        spit_cluster_m <- read.table(controlled_file_path, sep = "\t", header = TRUE, row.names = 1)
    } else if (file.exists(non_controlled_file_path)){
        spit_cluster_m <- read.table(non_controlled_file_path, sep = "\t", header = TRUE, row.names = 1)
    } else {
        print("No SPIT cluster matrices found. Please run SPIT dtu before running the cluster module.", flush=TRUE)
        return(NULL)
    }
    spit_cluster_m <- as.matrix(spit_cluster_m)
    spit_cluster_m <- ifelse(is.na(spit_cluster_m), 0, spit_cluster_m)
    if(!include_shared_dtu){
        rows_all_zeros <- which(rowSums(spit_cluster_m == 0) == ncol(spit_cluster_m))
        rows_all_ones <- which(rowSums(spit_cluster_m == 1) == ncol(spit_cluster_m))
        rows_to_extract <- union(rows_all_zeros, rows_all_ones)
        spit_cluster_m <- spit_cluster_m[-rows_to_extract, ]
    }
    if(all(rowSums(spit_cluster_m == 1) == ncol(spit_cluster_m))){
        print("All DTU events are shared amongst case samples, cannot apply clustering.", flush=TRUE)
        return(NULL)
    }
    distance <- dist(t(spit_cluster_m), method = "binary")
    set.seed(123)
    jitter_amount <- 0.01
    distance <- distance + abs(jitter(distance, amount = jitter_amount))
    clusters <- hclust(distance, method = "complete")
    dendrogram <- as.dendrogram(clusters)
    gene_cor_dist <- dist(spit_cluster_m, method = "binary")
    pdf_h <- nrow(spit_cluster_m) * 0.1 + 5
    pdf_w <- 11
    color_vector <- rep("White", ncol(spit_cluster_m))
    if (cov!=FALSE) {
        metadata = read.table(pheno, header=TRUE, sep = '\t')
        id_order <- colnames(spit_cluster_m)
        metadata_reordered <- metadata[match(id_order, metadata$id), ]
        cov_v <- metadata_reordered[cov]
        unique_cov <- unique(metadata[[cov]])
        num_colors <- length(unique_cov)
        color_palette <- viridis(num_colors)
        color_mapping <- setNames(color_palette, unique_cov)
        color_vector <- color_mapping[metadata[[cov]]]
    }
    pdf(file="spit_dendrogram.pdf", height=pdf_h, width=pdf_w)
    heatmap.2(spit_cluster_m,
    Colv = dendrogram,
    cexRow = 0.7, cexCol = 0.7,
    offsetRow = 0.01, offsetCol=0.01,
    colsep = c(1:ncol(spit_cluster_m)),
    rowsep = c(1:nrow(spit_cluster_m)),
    distfun = function(x) gene_cor_dist,
    col = colorRampPalette(brewer.pal(3, col_pallete))(2),
    ColSideColors = color_vector,
    linecol = '#E41653',
    lwid=c(1, 8),
    density.info = 'none',
    key.xtickfun = function() {
        breaks <- parent.frame()$breaks
        return(list(
        at = parent.frame()$scale01(c(breaks[1], breaks[length(breaks)])),
        labels = c(as.character(breaks[1]), as.character(breaks[length(breaks)]))
        ))
    },
    key.par=list(mgp=c(0.5, 0.5, 0),
    mar=c(pdf_h-3, pdf_w/5, pdf_h-3, pdf_w/5)),
    main = "Dendrogram of Case Samples Based on DTU Events Detected by SPIT")
    dev.off()
}
