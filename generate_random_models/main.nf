#!/usr/bin/env nextflow

params.template = "/media/houcem/theDrum/BILIM/nhlamulo/CYP2B6/CYP2B6_WT.pdb"
params.n = 100


process pepSeqStat {
	// requires EMBOSS package 

	output: 
		file("pepstats.dat") into pepstat

	"""
	pdb2fasta ${params.template} > template.fa 
	pepstats -sequence template.fa -outfile pepstats.dat 
	"""

}


process generateRandomSeq {
	echo true
	input: 
		file pep_stat_file from pepstat


	"""
	res_number=`grep "Residues =" $pep_stat_file |awk {'print \$7'}`
	echo \$res_number
	makeprotseq -length \$res_number \
		-outseq allseq.fa \
		-amount ${params.n} \
		-pepstatsfile pepstats.dat \
		-auto
	"""
}

process generatePirTemplate {
	echo true


		"""
		import modeller as *

		"""

}