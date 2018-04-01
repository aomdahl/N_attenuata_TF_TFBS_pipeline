# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 18:24:16 2018

@author: aomda
"""
import os
import subprocess as cmd
import argparse
import numpy as np



def CreateNewDir(par_path, name):
    if not os.path.exists(par_path +"/" + name):
        os.makedirs(par_path +"/" + name)    

#Gets the reference directory relative to aTf we want.
def GetReferenceDirectory(ref_dir, tf):
    dir_list = os.listdir(ref_dir)
    for r in dir_list:
        if tf in r:
            return ref_dir + "/" + r, r

            #Note that we need to keep track of the length of each term as well.

 #observe also the problem here... if there are multiple motifs in one file, this fails tragically.
"""def ReadInPWM(file_handle):
    with open(file_handle) as f_in:
        counter = 0
        all_weights = []
        for line in f_in:
            if ">" in line:
                counter += 1
                continue
            else:
                curr_weights = line.strip().split('\t')
                all_weights.append(curr_weights)
    return all_weights"""

def ParseMotifFile(file_handle):
    pwm_list = []
    with open(file_handle) as f_in:
        counter  = 0
        line_list = None
        for line in f_in:
            if ">" in line:
                if counter != 0:
                    pwm_list.append(line_list)
                counter += 1
                line_list = []
            else:
                curr_weights = line.strip().split('\t')
                line_list.append(curr_weights)
        pwm_list.append(line_list)
    return pwm_list
            
                



def AdjustPWMLengthsAll(pwms, min_len):
    for p in pwms:
        AdjustPWMLength(p, min_len)

def AdjustPWMLength(pwm, min_len):
    EQUAL = 0.25
    if len(pwm) < min_len:
        while len(pwm) < min_len:
            pwm.append([EQUAL,EQUAL,EQUAL,EQUAL])  #default all the values to equal probability so it will have no effect in averaging.
    elif len(pwm) > min_len:
        print "Unexpected condition, please revisit this case..."
    else:
        return

def AveragePWMs(pwm_list):
    
    print pwm_list
    raw_input()
    return np.average(np.array(pwm_list).astype(np.float), axis = 0)


if __name__ == '__main__':
        #Argument parsing
        parser = argparse.ArgumentParser(description='A function to average PWMs as they are passed in')
        parser.add_argument("pwm_paths", help = "A file containing all the paths to the single PWM files we wish to merge")
        args = parser.parse_args()

        pwm_super_list = []
        max_length = 0
        with open(args.pwm_paths) as p_in:
            for line in p_in:
               add_in = ParseMotifFile(line.strip())
               pwm_super_list += add_in
               if len(add_in) > max_length:
                   max_length = len(add_in)

        AdjustPWMLengthsAll(pwm_super_list, max_length)

        averaged_pwm  = AveragePWMs(pwm_super_list)
        print averaged_pwm