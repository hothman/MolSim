#!/usr/bin/python3

import argparse
import sys
import numpy as np
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
        return self.covarience_file.split(".cov")[0]

    def output(self, output_file, tag):
        output_list = [tag,self.S_vib]
        with open( output_file , 'w') as file:
            file.writelines( ','.join( output_list  )+"\n" )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=" The script extracts the vibrational entropy from ENCoM generated .cov file")
    parser.add_argument("--cov", help="covarience file")
    parser.add_argument("--output", help="output file (CSV)")
    parser.add_argument("--tag", help="tag to include in the output csv file")
    args = parser.parse_args()
    covarience = Encom(args.cov)
    covarience.getSvib()
    if args.tag == None: 
        tag = covarience.getTag()
    else: 
        tag = args.get
    covarience.output(args.output, tag)


