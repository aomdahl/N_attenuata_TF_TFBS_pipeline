#!/usr/bin/env bash
 
#This is the workflow used for both the RNA-seq and MA datasets.
#Here I only show some of the steps for convenience
#Extract the relevant files and condense them to the relevant directories using internal_motif_analysis.py
#This version, made on 3/29, was run on the unfiltered data (not chosen by conservation or anything)
mkdir nofilter_MA-MA_comparison_results
python3 internal_motif_analysis.py /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/ma.12.29.TF.superTable.remade.tsv /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/MA_results/ /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_MA-MA_comparison_results/

mkdir nofilter_MA-RNA_comparison_results
python3 internal_motif_analysis.py /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/ma.12.29.TF.superTable.remade.tsv /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNASeq_results/ /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_MA-RNA_comparison_results/

mkdir nofilter_RNA-RNA_comparison_results
python3 internal_motif_analysis.py /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNA.8.17.TF.superTable.tsv /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNASeq_results/ /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_RNA-RNA_comparison_results/

mkdir nofilter_RNA-MA_comparison_results
python3 internal_motif_analysis.py /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/RNA.8.17.TF.superTable.tsv /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/MA_results/ /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_RNA-MA_comparison_results/





#Check for any error in the files:
grep -r "--" ./no_filter_MA-MA_comparison_results/*
grep -r "--" ./no_filter_RNA-RNA_comparison_results/*
grep -r "--" ./no_filter_MA-RNA_comparison_results/*
grep -r "--" ./no_filter_RNA-MA_comparison_results/*

#Run the final script
bash motif_comparison.sh /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_RNA-RNA_comparison_results/
bash motif_comparison.sh /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_MA-RNA_comparison_results/
bash motif_comparison.sh /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_MA-MA_comparison_results/
bash motif_comparison.sh /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_RNA-MA_comparison_results/


#Now, we want to condense the results into a readable table. Here I use the Homer html parser I built over the summer:
mkdir nofilter_results
python3 HomerTools.py --dirs /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_MA-MA_comparison_results/ > /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_results/MA-MA-comparison_summary.tsv &
python3 HomerTools.py --dirs /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_RNA-RNA_comparison_results/ > /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_results/RNA-RNA-comparison_summary.tsv &
python3 HomerTools.py --dirs /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_MA-RNA_comparison_results/ > /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_results/MA-RNA-comparison_summary.tsv &
python3 HomerTools.py --dirs /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_RNA-MA_comparison_results/ > /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/nofilter_results/RNA-MA-comparison_summary.tsv &

