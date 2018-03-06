#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon March 1 18:24:16 2018

@author: aomda
"""

#$1 is the folder you want to step through


for f in $1/*; do
        if [[ -d $f ]]; then
                #do your thing. We want to only keep 
                cd $f
                compareMotifs.pl ./*filtered_for_comparison.motif ./*_ref.motif -reduceThresh 0.9 -pvalue 0.001
                cd ../            

fi
