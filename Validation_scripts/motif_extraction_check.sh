#!/usr/bin/env bash
#This tool will run the motif comparison of motifs against each other
#$1 is the folder you want to step through
cd $1
for d in ./* ; do
	if [[ -d $d ]]; then
		cd $d
		if [ -d "homerResults" ]; then
			cd ../
			continue
		fi
		if [ ! -s ./*filtered_for_comparison.motif ]; then
			cd ../
			continue
		fi
		if [ ! -s ./*_ref.motif ]; then
			cd ../
			continue
		fi	
		grep "--" ./*
		cd ../
	fi
done
