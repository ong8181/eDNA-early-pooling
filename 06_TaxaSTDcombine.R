####
#### No.6: STD sequences assignment
#### 2021.11.29 Ushio
#### R 4.1.2
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse") #1.3.1, 2021.11.11

# Generate output directory
od <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od, end = -3), "Out")); rm(od)
dir.create(output_folder)

# Load Claident taxa assignment
claident_tax <- read.delim("04_TaxaAssignmentOut/merge_classigntax")
head(claident_tax)

# Load Blastn taxa assignment for STD sequences
blastn_std0 <- read.table("05_identSTD_BLASTnOut/STDseqOut.txt")
#V1   qseqid      query or source (e.g., gene) sequence id
#V2   sseqid      subject  or target (e.g., reference genome) sequence id
#V3.  pident      percentage of identical matches
#V4.  length      alignment length (sequence overlap)
#V5.  mismatch    number of mismatches
#V6.  gapopen     number of gap openings
#V7.  qstart      start of alignment in query
#V8.  qend        end of alignment in query
#V9.  sstart      start of alignment in subject
#V10. send        end of alignment in subject
#V11. evalue      expect value
#V12. bitscore    bit score
cond1 <- blastn_std0$V5 < 5 & blastn_std0$V4 > 160 & blastn_std0$V6 < 1 # Alignment length and No. of mismatch and gap
blastn_std <- blastn_std0[cond1,]

# Check claident taxa assignments of the potential std sequences
potential_std_id <- match(blastn_std$V1, claident_tax$query)
claident_tax[potential_std_id, "family"]
claident_tax[potential_std_id, "species"]

# Replace STD taxa names with claident taxa assigment
claident_tax$species <- as.character(claident_tax$species)
claident_tax[potential_std_id, "species"] <- as.character(blastn_std$V2)
start_ncol <- which("order" == colnames(claident_tax))
end_ncol <- which("genus" == colnames(claident_tax))
claident_tax[potential_std_id, start_ncol:end_ncol] <- "STDseqs"

# Output new claident tax table
write.csv(claident_tax, sprintf("%s/claident_tax_revise.csv", output_folder), row.names = F)

# save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

