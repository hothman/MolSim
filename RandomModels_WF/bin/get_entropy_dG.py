#!/usr/bin/python3

import argparse
import sys
import re
import warnings
import subprocess
import os


class Encom:
    def __init__(self, covarience_file):
        self.covarience_file = covarience_file
        self.eigenvalues = []
        self.eigenvectors =  []

    def getSvib(self):
        with open(self.covarience_file, 'r') as file_to_read: 
            self.S_vib = file_to_read.readlines()[1].strip().split(":")[1]

    def getTag(self):
        self.tag = self.covarience_file.split(".cov")[0]

    def arrangeOutput(self):
        output_list = [self.tag,self.S_vib]
        return output_list


class Foldx:
    def __init__(self, foldx_output):
        self.foldx_output = foldx_output

    def getDg(self): 
        with open(self.foldx_output, 'r') as f: 
            lines = f.readlines()
        Dg = lines[0].split("\t")[1]
        internal_tag = lines[0].split("\t")[0].replace(".pdb", "").replace("./","")
        return [internal_tag,Dg]

def output(output_file, output_list, dg):
    s_vib_dg = output_list+dg
    with open( output_file , 'w') as file:
        file.writelines( ','.join( s_vib_dg  )+"\n" )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=" The script extracts the vibrational entropy from ENCoM generated .cov file")
    parser.add_argument("--cov", help="covarience file")
    parser.add_argument("--output", help="output file (CSV)")
    parser.add_argument("--tag", help="tag to include in the output csv file")
    parser.add_argument("--foldx", help="foldx fxout file" )
    args = parser.parse_args()
    covarience = Encom(args.cov)
    covarience.getSvib()
    if args.tag == None: 
        tag = covarience.getTag()
    else: 
        tag = args.get
    output_encom = covarience.arrangeOutput()
    foldx_processing = Foldx(args.foldx)
    output_foldx = foldx_processing.getDg()
    output(args.output, output_encom,output_foldx)


