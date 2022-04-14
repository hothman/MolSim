#!/usr/bin/python3

import argparse
import sys
import numpy as np
import re
import warnings
import subprocess
import os

path="/path/to/folder"

for file in os.listdir(os.chdir(path)):
    if file.endswith(".eigen"):
        class Encom:
            def __init__(self, eigenfile):
                self.eigenfile = eigenfile
                self.eigenvalues = []
                self.eigenvectors =  []

            def parse_eigen(self):
                with open(self.eigenfile, 'r') as file:
                    lines = file.read()
                raw_format = lines.split('\n\n')
                eigenvalues = [] ; eigenvectors = []
                for mode in raw_format:
                    if "Eigenvalue:" in mode:
                        header_and_data = re.compile('.*Eigenvalue:|\n').split(mode )
                        clean_list = [elem for elem in header_and_data if elem]
                        eigenvalues.append( float(clean_list[0]) )
                        current_eigenvector = []
                        for eigenvector_line in header_and_data[1:]:
                            component = eigenvector_line.split('\t')
                            try :
                                current_eigenvector.append( float(component[2]) )
                            except:
                                pass
                        eigenvectors.append( current_eigenvector )
                self.eigenvalues = eigenvalues
                self.eigenvectors = eigenvectors
            def calRMSIP(subspace1, subspace2, n = 100):
                if (isinstance(subspace1, Encom) and isinstance(subspace2, Encom)):
                    subspace1.parse_eigen()
                    sub1 = subspace1.eigenvectors
                    subspace2.parse_eigen()
                    sub2 = subspace2.eigenvectors
                    sum_outer = 0
                    for eigenvector1 in sub1[6:6+n]:
                        for eigenvector2 in sub2[6:6+n]:
                            sum_outer = sum_outer + (np.dot(np.array( eigenvector1 ), np.array( eigenvector2)) )**2
                    return np.sqrt(np.array(sum_outer) /float(n) )
                else:
                    raise TypeError ('Arguments must be instances of class Encom')    

            def getEntropy(wt_subspace, mut_subspace):
                if (isinstance(wt_subspace, Encom) and isinstance(mut_subspace, Encom)):
                    wt_subspace.parse_eigen()
                    eigen1 = wt_subspace.eigenvalues
                    mut_subspace.parse_eigen()
                    eigen2 = mut_subspace.eigenvalues
                    assert len(eigen2) == len(eigen1)
                    fraction_eigenvalues = np.array(eigen2[6:])/np.array(eigen1[6:])
                    return np.log(np.prod(fraction_eigenvalues))

            mut_encom = Encom(file)
            mut_encom.parse_eigen()
            mut_eigenvectors = mut_encom.eigenvectors
            mut_eigenvalues = mut_encom.eigenvalues
            wt_encom = Encom("/path/to/WT/.eigen")
            wt_encom.parse_eigen()
            mut_eigenvectors = wt_encom.eigenvectors
            mut_eigenvalues = wt_encom.eigenvalues
            rmsip= str(calRMSIP(mut_encom, wt_encom))
            entropy = str( round(getEntropy(mut_encom, wt_encom), 3))

            print("entropy =", entropy, "kcal/mol", "rmsip =", rmsip)
