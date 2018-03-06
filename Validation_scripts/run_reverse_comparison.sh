#!/usr/bin/env bash
 
#This is the workflow used for both the RNA-seq and MA datasets.
#Here I only show some of the steps for convenience
#Extract the relevant files and condense them to the relevant directories using internal_motif_analysis.py

python3 internal_motif_analysis.py /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/TopRNATFCandidates_1_6_18.tsv /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNASeq_results/ /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNA-RNA_comparison_results/

#Check for any error in the files:
grep -r "--" ./RNA-RNA_comparison_results/*

#Run the final script
bash motif_comparison.sh /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNA-RNA_comparison_results/

#we are also interested, for curiosity's sake, in running this in reverse on RNA ids to MA networks, and visa versa
python3 internal_motif_analysis.py /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/TopRNATFCandidates_1_6_18.tsv /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/MA_results/ /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNA-MA_comparison_results/

python3 internal_motif_analysis.py /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/TopMATFCandidates_1_6_18.tsv /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNASeq_results/ /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/MA-RNA_comparison_results/

bash motif_comparison.sh /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNA-MA_comparison_results/
bash motif_comparison.sh /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/MA-RNA_comparison_results/

#Now, we want to condense the results into a readable table. Here I use the Homer html parser I built over the summer:
python3 HomerTools.py --dirs /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/MA-MA_comparison_results/ > /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/MA-MA-comparison_summary.tsv &
python3 HomerTools.py --dirs /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNA-RNA_comparison_results/ > /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNA-RNA-comparison_summary.tsv &
python3 HomerTools.py --dirs /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/MA-RNA_comparison_results/ > /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/MA-RNA-comparison_summary.tsv &
python3 HomerTools.py --dirs /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNA-MA_comparison_results/ > /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNA-MA-comparison_summary.tsv &

