#!usr/bin/env nextflow 

ligands = Channel.fromPath("${params.mol2}/*.mol2")
myreceptor = Channel.fromPath("${params.myreceptor}/*.pdb")
myreceptor2 = Channel.fromPath("${params.myreceptor}/*.pdb")
myreceptor2.subscribe { println it }
config_file = Channel.fromPath("${params.PLANTS_config}")
ligands2 = Channel.fromPath("${params.mol2}/*.mol2")
//ligands.subscribe { println it } 

process prepare_conformer {
	input:
		file receptor from myreceptor
	output:
		file("*.mol2") into receptor_mol 
	script: 
		name = receptor.baseName.replaceFirst(".pdb","")

	"""
	SPORES_64bit --mode complete $receptor ${name}.mol2 
	rm *_bad.mol2
	"""
}

process docking {
	publishDir './CLC3VS' , mode: 'move', overwrite: true
	maxForks 10
	input: 
		each file(ligand) from ligands
		each file(receptor) from receptor_mol
		file config from config_file
	output: 
		file "*mol2VSdir*" into VSfolder 

	"""
	SPORES_64bit --mode protstates $ligand ligandprotonated.mol2
	sed -i 's/XXXXX/$receptor/' $config
	sed -i 's/YYYYY/ligandprotonated.mol2/' $config
	sed -i 's/ZZZZZ/${receptor}_${ligand}VSdir/' $config
	PLANTS1.2_64bit --mode screen $config
	"""

}

