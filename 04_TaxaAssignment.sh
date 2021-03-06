####
#### Commands to demultiplex bcl to fastq
####

#----------- Taxa assignment using claident -----------#
FASTA_FILE="OTU_seqs.fa"
SEQ_FOLDER="03_OTUClusteringOut"
OUTPUT_FOLDER="04_TaxaAssignmentOut"
mkdir ${OUTPUT_FOLDER}
cd ${SEQ_FOLDER}

# Check overall_class
clmakecachedb --blastdb=overall_class --numthreads=36 ${FASTA_FILE} overall_class_cache
clidentseq --blastdb=overall_class_cache --numthreads=36 ${FASTA_FILE} overall_class_clidentseq
classigntax --taxdb=overall_class --maxpopposer=0.05 --minsoratio=19 overall_class_clidentseq overall_class_classigntax

# Overall_genus
clmakecachedb --blastdb=overall_genus --numthreads=36 ${FASTA_FILE} overall_genus_cache
clidentseq --blastdb=overall_genus_cache --numthreads=36 ${FASTA_FILE} overall_genus_clidentseq
classigntax --taxdb=overall_genus --maxpopposer=0.05 --minsoratio=19 overall_genus_clidentseq overall_genus_classigntax

# Merge identification results (overall_class + overall_genus)
clmergeassign --priority=equal --preferlower overall_genus_classigntax overall_class_classigntax merge_classigntax

# Delete large files
rm -r overall_class_cache
rm -r overall_genus_cache

# Move file
mv overall_class_clidentseq ../${OUTPUT_FOLDER}
mv overall_class_classigntax ../${OUTPUT_FOLDER}
mv overall_genus_clidentseq ../${OUTPUT_FOLDER}
mv overall_genus_classigntax ../${OUTPUT_FOLDER}
mv merge_classigntax ../${OUTPUT_FOLDER}
