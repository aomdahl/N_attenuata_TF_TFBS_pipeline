# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 18:24:16 2018

@author: aomda
"""
import os
import subprocess as cmd
import argparse

#Pass_in is a line from the extracted file
#returns the identfier, as well as the binding motif and the associated TF
def extractDirIDKey(pass_in):
    core_gene = pass_in[11]
    cg_fam = pass_in[12]
    return cg_fam + "-" + core_gene, pass_in[4]

def ExtractTFKey(pass_in):
    tf_name = pass_in[0].replace(".t1","")
    return tf_name

#Extracts the wanted motif information from a file and puts it in a mergable file elswhere
def combineMotifFilesBySeq(motif_seq, grep_file, outfile):
    extract_length = str(len(motif_seq))
    args_list =["grep", "-A", extract_length, motif_seq, grep_file, ">>", outfile]

    try:
        os.system(" ".join(args_list))
        #add a newline to the end if it doesn't have one already
        os.system("sed -i -e \'a\\\' " + outfile)
    except:
        print("unxpected error on" + grep_file)

def combineMotifFiles(thing_one, thing_two, outfile):
    try:
        os.system("cat " +  thing_one +">> " + outfile)
        os.system("cat " +  thing_two +">> " + outfile)
    except:
        print("unexpected error on" + outfile)


#Tracking set is used to see which queries have already been made.
# This updates the list of queries already made  
def updateTrackingSet(search_key, motif, tf, ts):
    ts.add(search_key+ "_" + motif + "_" + tf)
    

#this checks thelist of queries already made, and returns true or false
def InTrackingSet(search_key, motif, tf, ts):
    if search_key + "_" + motif + "_" + tf in ts:
        return True
    return False

def CreateNewDir(par_path, name):
    if not os.path.exists(par_path +"/" + name):
        os.makedirs(par_path +"/" + name)    

#Gets the reference directory relative to aTf we want.
def GetReferenceDirectory(ref_dir, tf):
    dir_list = os.listdir(ref_dir)
    for r in dir_list:
        if tf in r:
            return ref_dir + "/" + r, r

        


if __name__ == '__main__':
        #Argument parsing
        parser = argparse.ArgumentParser(description='A quick function to check for similar motifs based on TFs')
        parser.add_argument("dataInput", help="Specify the path to the filtered TFs")
        parser.add_argument("dataDir", help = "Specify the directory with the found TF data")
        parser.add_argument("outFileDest", help= " The destination directory for output files")
        parser.add_argument("-r", "--refDir",help = "Specify the reference directory to compare the found motifs to")
        args = parser.parse_args()

        if args.refDir is None:
            args.refDir = args.dataDir


        ##Set the search database
        tracking_set = set()
        dir_list = os.listdir(args.dataDir)
        with open(args.dataInput) as f_in:
            for line in f_in:
                line = line.strip()
                if "PutativeTF" not in line:
                    search_key, p_motif = extractDirIDKey(line.split('\t'))
                    p_tf = ExtractTFKey(line.split('\t'))
                    dir_name_out = search_key + "_against_" + p_tf
                    #If we are just going to extract the same data again...
                    if InTrackingSet(search_key, p_motif, p_tf, tracking_set):
                        #skip the step
                        continue

                    if search_key in dir_list: 
                        CreateNewDir(args.outFileDest, dir_name_out)
                        #grep out the motifs that we want, and put it in a file
                        outfile_name = args.outFileDest + "/" +dir_name_out +"/"
                        outfile_name_query = outfile_name + search_key + "_filtered_for_comparison.motif"
                        file_name_1kb = args.dataDir +"/" + search_key + "/" + "1kb_topMatches.motif"
                        file_name_2kb = args.dataDir +"/" + search_key + "/" + "2kb_topMatches.motif"
                        combineMotifFilesBySeq(p_motif, file_name_1kb, outfile_name_query)
                        combineMotifFilesBySeq(p_motif, file_name_2kb, outfile_name_query)
                        #now get the information for the reference motifs
                        try:
                            ref_dir, tf_id = GetReferenceDirectory(args.refDir, p_tf)
                        except TypeError:
                            #this means that the modele was not made
                            print ("No module data exists for", p_tf)
                            continue

                        outfile_name_ref = outfile_name + "/" + tf_id + "_ref.motif"
                        file_name_1kb_ref = ref_dir +"/" + "1kb_topMatches.motif"
                        file_name_2kb_ref = ref_dir +"/" + "2kb_topMatches.motif"

                        combineMotifFiles(file_name_1kb_ref,file_name_2kb_ref, outfile_name_ref)

                        
                    else:
                        print (search_key + " dir doesn't exist. Moving on.")
                    updateTrackingSet(search_key, p_motif, p_tf, tracking_set)               
        
#This step will run on its own-- too much to try and do it all here. Let's just get it primed here.
#compareMotifs.pl /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/MA-MA_comparison_results/AP2-EREBP-NIATv7_g01614_filtered_for_comparison.motif /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/MA-MA_comparison_results/ -known /mnt/d/AshtonData/MPI_all_data/SummaryFiles_8-17/MA_results/HB-NIATv7_g14998/2kb_topMatches_FORTESTING.motif -reduceThresh 0.8 -pvalue 0.001