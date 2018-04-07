#!/usr/bin/env bash
#This tool will run the motif comparison of motifs against each other
#$1 is the folder you want to step through
cd $1
for d in ./* ; do
	if [[ -d $d ]]; then
		cd $d
		if [ ! -s ./*filtered_for_comparison.motif ]; then
			cd ../
			rm -r $d
			continue
		fi
		if [ ! -s ./*_ref.motif ]; then
			cd ../
			rm -r $d
			continue
		fi	
		echo $d
		sed -i '/^$/d' ./*filtered_for_comparison.motif
		sed -i '/^$/d' ./*_ref.motif
		compareMotifs.pl ./*filtered_for_comparison.motif ./ -known ./*_ref.motif -reduceThresh 0.9
		cd ../
	fi
done
