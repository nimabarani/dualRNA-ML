library(DESeq2)
library(logr)

set.seed(42)

get_count_table <- function(path, experiment, pattern) {
    i <- 1
    sample <- 'uninfected'
    if (experiment == 'infection')
        sample <- 'infected'
    
    col_cases <- c(NA, 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', NA)
    sra_dir <- file.path(path, experiment)
    files_list <- list.files(
        path = file.path(path, experiment),
        pattern = pattern
    )
    
    df <- lapply(file.path(sra_dir, files_list), function(file) {
        tmp_df <- read.csv(file, sep = '\t', colClasses = col_cases, row.names = 1, skip = 1)
        colnames(tmp_df) <- paste0(sample, '_rep_', i)
        i <<- i+1
        return(tmp_df)
    })
    
    return(df)
}

get_results <- function(treated_df, untreated_df) {
    read_counts <- data.frame(untreated_df, treated_df)
    meta <- data.frame(
        condition = c(
            rep('untreated', length(untreated_df)),
            rep('treated', length(treated_df))
        ),
        row.names = colnames(read_counts)
    )
    ## Create DESeq2Dataset object
    dds <- DESeqDataSetFromMatrix(
        countData = read_counts,
        colData = meta,
        design = ~condition
    )
    dds$condition <- relevel(dds$condition, ref = 'untreated')
    dds <- DESeq(dds)
    # Get log2 fold change, base mean, and p-value for each gene
    res <- results(dds)
    # Get normalized counts for each sample and gene
    norm_counts <- counts(dds, normalized=TRUE)
    # Combine normalized counts, log2 fold change, base mean, and p-value into a single data frame
    results_df <- as.data.frame(cbind(norm_counts, res$log2FoldChange, res$baseMean, res$pvalue, res$padj))
    # Rename columns
    colnames(results_df) <- c(colnames(norm_counts), "log2FoldChange", "baseMean", "pvalue", 'padj')
    # Log: Shape of results
    log_print('Shape of results:', hide_notes = T, blank_after = F)
    log_print(dim(results_df))
    return(results_df)
}

get_filtered_df <- function(results) {
    notnull_res <- results[complete.cases(results), ]
    notnull_res_copy <- data.frame(notnull_res)
    # Log: Shape of non null results
    log_print('Shape of non nulls:', hide_notes = TRUE, blank_after = FALSE)
    log_print(dim(notnull_res), hide_notes = TRUE)
    
    if (nrow(notnull_res[notnull_res_copy$padj < 0.1 & notnull_res_copy$log2FoldChange >= 1, ])) {
        notnull_res[notnull_res_copy$padj < 0.1 & notnull_res_copy$log2FoldChange >= 1, ]$log2FoldChange <- 'up'
    }
    
    if (nrow(notnull_res[notnull_res_copy$padj < 0.1 & notnull_res_copy$log2FoldChange <= -1, ])) {
        notnull_res[notnull_res_copy$padj < 0.1 & notnull_res_copy$log2FoldChange <= -1, ]$log2FoldChange <- 'down'
    }
    
    if (nrow(notnull_res[notnull_res_copy$padj > 0.1, ])) {
        notnull_res[notnull_res_copy$padj > 0.1, ]$log2FoldChange <- 'nd'
    }
    
    filtered_df <- notnull_res[notnull_res$log2FoldChange %in% c('up', 'down', 'nd'), ]
    
    # Log: Shape of filtered results
    log_print('Shape of filtered results:', hide_notes = T, blank_after = F)
    log_print(dim(filtered_df), hide_notes = T)
    
    # Log: Number of unique values
    log_print('Number of unique values:', hide_notes = T, blank_after = F)
    log_print(table(filtered_df$log2FoldChange), hide_notes = T)
    return(filtered_df)
}

# ---------------------------------------------------------
# Set the directory where the CSV files are stored
parent_dir <- '~/Developer/counts'

# Get a list of all subdirectories in the parent directory
sub_dirs <- list.dirs(parent_dir, recursive = FALSE)

for (dir_path in sub_dirs) {
    tmp <- file.path(dir_path, 'main.log')
    lf <- log_open(tmp, show_notes = F)

    log_print(basename(dir_path), hide_notes = T, blank_after = F)

    # Host
    # Log: Host
    log_print('Host:', hide_notes = T, blank_after = F)
    untreated_hosts_df <- get_count_table(path = dir_path, experiment = 'host', pattern = '\\_host.txt$')
    treated_hosts_df <- get_count_table(path = dir_path, experiment = 'infection', pattern = '\\_host.txt$')
    host_res <- get_results(
        treated_df = treated_hosts_df,
        untreated_df = untreated_hosts_df
    )
    host_filtered_res <- get_filtered_df(host_res)
    write.csv(host_filtered_res, file.path(dir_path, 'host.csv'))

    # Log: Pathogen
    log_print('Pathogen:', hide_notes = T, blank_after = F)
    untreated_pathogen_df <- get_count_table(path = dir_path, experiment = 'pathogen', pattern = '\\_pathogen.txt$')
    treated_pathogen_df <- get_count_table(path = dir_path, experiment = 'infection', pattern = '\\_pathogen.txt$')
    pathogen_res <- get_results(
        treated_df = treated_pathogen_df,
        untreated_df = untreated_pathogen_df
    )
    pathogen_filtered_res <- get_filtered_df(pathogen_res)
    write.csv(pathogen_filtered_res, file.path(dir_path, 'pathogen.csv'))

    log_close()
}
