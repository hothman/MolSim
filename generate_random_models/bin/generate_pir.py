#!/usr/bin/env python

from Bio import SeqIO
import argparse
from modeller import *

def getPirTarget(fasta_seq_file):
	"""
	read fasta file and generate a PIR format of 
	unaligned sequences by MODELLER
	"""
	env = Environ()
	env.libs.topology.read('${LIB}/top_heav.lib')
	env.libs.parameters.read('${LIB}/par.lib')
	aln = Alignment(env)
	# append sequences to alm object
	for i, record in enumerate(SeqIO.parse(fasta_seq_file, "fasta") ):
		print(i)
		aln = Alignment(env)
		aln.append_sequence(str( record.seq.upper()) )
		aln[0].code = "RANDOMSEQ"+str(i)
		aln.write(file="./RANDOMSEQ"+str(i)+".pir")


def getPirTemplate(pdb_file, output_file):
	"""
	generates PIR file for a PDB file
	"""
	env = Environ()
	mdl = Model(env, file=pdb_file)
	for c in mdl.chains:
		if c.filter( structure_types='structureN structureX structureM', chop_nonstd_termini=True):
			print(c)
			atom_file, align_code = c.atom_file_and_code(output_file)
			c.write(output_file, atom_file, align_code, format='PIR', chop_nonstd_termini=True)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=" Generate PIR files for MODELLER")
	parser.add_argument("--mode", help="template or target")
	parser.add_argument("--fasta", help="input fasta file")
	parser.add_argument("--output", help=" Output name of PIR file")
	parser.add_argument("--pdb", help=" PDB file of template")
	args = parser.parse_args()

	if args.mode == "target":
		getPirTarget(args.fasta)
	elif args.mode == "template":
		getPirTemplate(args.pdb, args.output)
