#!/usr/bin/env sh
#Takes the path to the module file and a name to give the module
#$1 is the module file path
#$2 is the module name
alias compariMotif='python /home/likewise-open/ICE/aomdahl/bin/slimsuite/tools/comparimotif_V3.py'
alias BLASTpSearch='python /home/likewise-open/ICE/aomdahl/Toolkit/MotifRecognition/BLASTPSearch.py'
alias HOMER2PRESTO='python /home/likewise-open/ICE/aomdahl/Toolkit/MotifRecognition/HOMER2PRESTO.py'
alias motifDBSearch='python /home/likewise-open/ICE/aomdahl/Toolkit/MotifRecognition/motifDBSearch.py'
alias tabularResults='python /home/likewise-open/ICE/aomdahl/Toolkit/MotifRecognition/tabularResults.py'
alias buildHOMERScript='python /home/likewise-open/ICE/aomdahl/Toolkit/MotifRecognition/build_HOMER_script.py'
alias getModuleOrthologs='python /home/likewise-open/ICE/aomdahl/Toolkit/MotifRecognition/getModuleOrthologs.py'
alias topMotifs='python /home/likewise-open/ICE/aomdahl/Toolkit/MotifRecognition/top_motif.py'

#Get the promoter sequences for the input file in the folder

mkdir $2
cp $1 $2
cd $2

mkdir 1kb
mkdir 2kb

##MAKE IT STOP IF WE DON't get valid input....
#Extract the promoter sequences we want
perl ~/Toolkit/Get_seq_fasta.pl -i *_module.txt -o ./1kb/1kb_promoters.fa -d ~/Datasets/GenomeData/1kb_promoter_sequences.fa -c 1
perl ~/Toolkit/Get_seq_fasta.pl -i *_module.txt -o ./2kb/2kb_promoters.fa -d ~/Datasets/GenomeData/2kb_promoter_sequences.fa -c 1
read -p $'Did the output come out valid??'

###References
POSITIVES=_promoters.fa
PROMOTER_PATH=~/Datasets/GenomeData/
NEGATIVES=_promoter_sequences.fa
 

CISBPLIB=/home/likewise-open/ICE/aomdahl/Datasets/CIS-BP/knownMotifDatabase.motifs
PPLIB=/home/likewise-open/ICE/aomdahl/Datasets/PlantPan/knownMotifDatabase.motifs
TOMATO_ORTHOLOG=/home/likewise-open/ICE/aomdahl/Datasets/GenomeData/NIAT_Orthologs/NIAT-tomato-mapping.tsv
###################
#For the different kb lengths


for folder in */; do
	cd $folder
	#Construct the HOMER script to run on both the positive and bootstrap samples
	build_HOMER_script $folder$POSITIVES $PROMOTER_PATH$folder$NEGATIVES -p 2 -l 8
	build_HOMER_script $folder$POSITIVES $PROMOTER_PATH$folder$NEGATIVES -p 2 -l 10


	#Loop through the Motif Identification folders with putative motifs
	for file in MotifIdentification*/; do
		echo "Starting HOMER analysis for $file"
		cd $file
		chmod 777 $file/MotifIdentificationScript*  #access the script
		perl /home/likewise-open/ICE/aomdahl/Toolkit/Run_shells_parallel.pl -n 9 -s MotifIdentificationScript*.sh -p Y  #Make the call

		#Take the HOMER output and convert into into PRESTO format for motif comparison
		echo "Starting bootstrap analysis..."
		HOMER2PRESTO -b ./PutativeMotifs/BootstrapMotifResults/ -o ./bootstrapPromoters.presto -a bootstrap -c 2
		echo "Starting positive data analysis..."
		HOMER2PRESTO -i ./PutativeMotifs/experimentalResults/homerMotifs.* -o ./positivePromoters.presto -a pTFBS -c 1e-5

		#Motif comparison
		compariMotif motifs=positivePromoters.presto searchdb=bootstrapPromoters.presto dna=T unmatched=TRUE

		#Keep those motifs that are least similar to those in the random bootstrap samples (judged by similarity score)
		topMotifs *compare.tdt positivePromoters.presto bootstrapPromoters.presto *.unmatched.txt ./PutativeMotifs/experimentalResults/homerMotifs* -m $2 #also include the homer reference file.
	
		cd ../ #Return to the original directory
		read -rsp $'Finished run of motifs in $file...\n' -n1 key
	done
	
	##Identify which genes from the module in which motifs actually appear for comparison checking.
	scanMotifGenomeWide.pl topMatches.motif *kb_promoters.fa -bed > topMotifSites.bed 

	#Split the promoter sequence files up within the file
	mkdir SplitPromoters
	perl  /home/likewise-open/ICE/aomdahl/Toolkit/split_multifasta.pl --input_file $folder$POSITIVES --output_dir ./SplitPromoters/

	#Conservation testing step-- using the default of 1000 bootstrap runs.
	getModuleOrthologs SplitPromoters topMotifSites.bed /home/likewise-open/ICE/aomdahl/Datasets/GenomeData/NIAT_Orthologs/SplitPromoterSequences/$folder/ $TOMATO_ORTHOLOG
	
	read -rsp $'Check the ortholog comparison output...\n' -n1 key
	
	#Starting reference TF lookup
	mkdir ReferenceTFs
	mkdir ReferenceTFs/PP
	mkdir ReferenceTFs/CIS-BP
	
	#Looking up identified motifs against the library
	compareMotifs.pl ConservationTesting/conservedMotifs.motifs ./ReferenceTFs/CIS-BP -known $CISBPLIB -matchThresh 0.8 -reduceThresh 0.9 -cpu 15  #CISBP
	compareMotifs.pl ConservationTesting/conservedMotifs.motifs ./ReferenceTFs/PP -known $PPLIB -matchThresh 0.8 -reduceThresh 0.9 -cpu 15 #PlantPan

	#This compares our motifs to the known motif database.  Using default species
	motifDBSearch ./ReferenceTFs

	#Print an output table so far-- preBlast result
	#This call still isn't working....
	#tabularResults -o TFGeneID_matches_preBLAST.tsv -t PreBlast

	#Pair them up with the blastp results, conservation pvalue 1e-5
	BLASTpSearch -cp 1e-5
	tabularResults -o blastResults.tsv -t PostBlast  #Print out all the results

	#Print the result library files.
	tabularResults -o pTF_library.tsv -t pTFLib

	tabularResults.py -t moduleResultsLib
done

