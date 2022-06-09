library(org.Hs.eg.db)

qld_data <- readRDS("/Users/jess/Work/2021_endo/endo_molecular_model/data/tidy_data/qld_counts_list.rds")

# Load other data to save everything in one RData file
exprs <- readRDS("../geo_analysis/cache/geo_rna_exprs.rds")
sample_df <- readRDS("../geo_analysis/cache/geo_rna_sample_df.rds")

# ensembl <- qld_data$gene_info$ensembl_id
ensembl <- rownames(exprs)
length(ensembl)

cols <- c("SYMBOL", "ENTREZID", "GENENAME")
df <- select(org.Hs.eg.db, keys=ensembl, columns=cols, keytype="ENSEMBL")
nrow(df)

# If duplicate entries, get the first entry
dup <- duplicated(df$ENSEMBL)
table(dup)

# dup_ensembl <- df[dup, "ENSEMBL"]
# df %>% filter(ENSEMBL %in% dup_ensembl) %>% nrow
# View(df %>% filter(ENSEMBL %in% dup_ensembl))

df <- df[! dup,]

# stopifnot(qld_data$gene_info$ensembl_id == df$ENSEMBL)
# gene_info <- cbind(qld_data$gene_info, df)
stopifnot(ensembl == df$ENSEMBL)
gene_info <- cbind(ensembl_id=ensembl, df)

# org.Hs.eg.db symbols and entrez ids are more complete than the old ones from the qld counts data
gene_info <- gene_info[,c("ensembl_id", "ENTREZID", "SYMBOL", "GENENAME")]
colnames(gene_info) <- c("ensembl_id", "entrez_id", "symbol", "gene_name")

# Reorder gene_info
m <- match(gene_info$ensembl_id, rownames(exprs))
gene_info <- gene_info[m,]

# Save
save(gene_info, exprs, sample_df, file="data/data.RData", version=2)
