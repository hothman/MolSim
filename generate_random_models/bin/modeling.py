#!/usr/bin/env python

import argparse
from modeller import *
from modeller.automodel import * 

parser = argparse.ArgumentParser(description=" create simple model")
parser.add_argument("--ali", help="alignment file")
args = parser.parse_args()

env = Environ()

## read PIR alignment file
aln = Alignment(env)
input = modfile.File(args.ali, 'r')
codes = []
while aln.read_one(input, alignment_format='PIR'):
    codes.append(aln[0].code)

print(codes)

env.io.atom_files_directory = ['.']
a = AutoModel(env, alnfile  = args.ali , knowns= codes[0], sequence = codes[1])
a.starting_model= 1                 
a.ending_model  = 1                                               
a.make()  
                  