# Load functions
first_col_to_rownames <- function(df) {
  rownames(df) <- df[[1]]
  df <- df[, -1]
  return(df)
}

# Read in RNA values where columns are unique sample IDs and rows are unique gene IDs
# 1. Convert to data frame
# 2. Move first column to rownames
expr_data_raw <- data.table::fread(expr_dir) %>%
  as.data.frame() %>%
  first_col_to_rownames()

# Limit expression data based on qualifier variable:
# 1. Select all sample IDs that contain qualifier, then grab the whole dataframe
# 2. Convert the first column to rownames
expr_data_sel <- expr_data_raw %>%
  dplyr::select(1 | contains(qual))

# Preprocess RNA data:
# 1. Transpose data to move genes to column names
# 2. Convert data to a data frame
expr_data_prep <- expr_data_sel %>%
  t() %>%
  as.data.frame()

# Filter genes for low expression
datExpr <- expr_data_prep %>%
  cleanRNAseq::filter_low_genes(min_expr, percent_cutoff,
                                metric)

# Grab sample IDs to filter the behavioral dataframe
sample_ids <- datExpr %>%
  rownames()

## Read in behavioral data where columns are behaviors and rows are unique sample IDs
# 1. Convert data.table to data.frame
# 2. Filter the first column (sample IDs), based on sample IDs in datExpr
datTraits <- data.table::fread(behav_dir) %>%
  as.data.frame() %>%
  dplyr::filter(.[[1]] %in% sample_ids) %>%
  first_col_to_rownames()

## Make sure sample IDs for both dataframes are identical
if (identical(sort(rownames(datExpr)), sort(rownames(datTraits)))) {
  message("You're all set! Sample IDs are consistent between datExpr and datTraits.")
} else {
  warning("Sample IDs do NOT match between datExpr and datTraits. Please double check your original data to ensure the number of samples are the same and that the sample IDs are identical.")
}

## Clear out garbage to save memory
rm(expr_data_prep, expr_data_raw, expr_data_sel)