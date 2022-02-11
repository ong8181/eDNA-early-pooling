####
#### F02. Figure Helper functions
####

# taxa name summarize function
taxa_name_summarize <- function(ps_object, taxa_rank, top_taxa_n = 10, taxa_always_include = NULL){
  tax_df <- as.data.frame(tax_table(ps_object))
  if(is.null(tax_df$rep_tax)) tax_df$rep_tax <- "Undetermined"

  # Search Others and Undetermined taxa
  tax_col1 <- which(colnames(tax_df) == taxa_rank) # Target rank
  tax_col2 <- which(colnames(tax_df) == "species") # Highest resolution
  ## Search unidentified taxa (= taxa is null from the target rank to the highest resolution)
  rep_tax_cond1 <- tax_df[,taxa_rank] == "" & !is.na(tax_df[,taxa_rank])
  rep_tax_cond2 <- apply(tax_df[,tax_col1:tax_col2] == "", 1, sum) == (tax_col2 - tax_col1) + 1
  
  # Replace taxa names
  tax_df[!rep_tax_cond1, "rep_tax"] <- as.character(tax_df[!rep_tax_cond1, taxa_rank])
  tax_df[rep_tax_cond1 & !rep_tax_cond2, "rep_tax"] <- "Others"
  
  # Re-import phyloseq object with revised tax_table
  ps_object2 <- phyloseq(otu_table(ps_object), sample_data(ps_object), tax_table(as.matrix(tax_df)))
  
  # Replace low abundance taxa name with Others
  taxa_abundance_rank <- aggregate(taxa_sums(ps_object2), by = list(tax_table(ps_object2)[,"rep_tax"]), sum)
  taxa_abundance_rank <- taxa_abundance_rank[order(taxa_abundance_rank$x, decreasing = T),]
  taxa_top <- taxa_abundance_rank[1:top_taxa_n,]
  
  if(is.null(taxa_always_include)) {
    include_taxa <- as.character(taxa_top[,1])
  } else {
    include_taxa <- unique(c(as.character(taxa_top[,1]), taxa_always_include))
  }
  
  low_tax <- is.na(match(tax_table(ps_object2)[,"rep_tax"], include_taxa))
  tax_table(ps_object2)[low_tax,"rep_tax"] <- "Others"
  
  return(ps_object2)
}
