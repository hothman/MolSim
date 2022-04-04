#!/usr/bin/env nextflow

params.template = ""
params.n = 100
params.env=""  // use this virtual envirement that contains all the dependencies
params.outfolder="./output"

process pepSeqStat {
	// requires EMBOSS package 

	output: 
		file("pepstats.dat") into pepstat

	"""
	pdb2fasta ${params.template} > template.fa 
	pepstats -sequence template.fa -outfile pepstats.dat 
	"""

}

process getTemplatePIR {
	// this process creates a pir file from the PDB file
	conda params.env
	output: 
		file "template.pir" into template_PIR
		file "template*.pdb" into template_structure
	
	"""
	generate_pir.py --mode template --pdb ${params.template} \
		--output template.pir
	"""
}


process generateRandomSeq {
	input: 
		file pep_stat_file from pepstat
	output: 
		file "allseq.fa" into all_random_seq

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


process getTargetPIR {
	// converts the randomly generated fasta sequences 
	// to separate PIR files
	input: 
		file(fasta) from all_random_seq
	output:
		file("*.pir") into target_pir

	"""
	generate_pir.py --mode target --fasta $fasta 
	"""
}


process modelRandomSeq {
	conda params.env

	publishDir params.outfolder, mode:'copy'
	input:
		each file(target) from target_pir
		each file(template_pdb) from template_structure
		file(template) from template_PIR
	output: 
		file("*.pdb")

	script:
		name = target.baseName.replaceFirst(".pir","")
	"""
	cat  $template > ${name}.ali
	cat $target >> ${name}.ali
	modeling.py --ali ${name}.ali
	mv ${name}.*.pdb ${name}.pdb
	"""

}