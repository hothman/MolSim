#!usr/bin/env nextflow 

ligands = Channel.fromPath("${params.mol2}/*.mol2")
receptors = Channel.fromPath("${params.receptors}/*.pdb")

config_file = Channel.fromPath("${params.PLANTS_config}")

ligands.subscribe { println it } 

process prepare_conformer {
	input:
		file receptor from receptors

	output:
		file("*.mol2") into receptor_mol 

	script: 
		name = receptor.baseName.replaceFirst(".pdb","")

	"""
	SPORES_64bit --mode complete $receptor ${name}.mol2 
	"""
}


process docking {
	input: 
		file receptor from receptor_mol
		file ligand from ligands
		file config from config_file
	output:
		file "${receptor_name}_${ligand_name}" into couple_rec_lig

		publishDir "virtual_screening", mode: 'copy'

	script:
		receptor_name = receptor.baseName.replaceFirst(".mol2","")
		ligand_name = ligand.baseName.replaceFirst(".mol2","")

	"""
	sed -i 's/XXXXX/${receptor_name}.mol2/' $config  
	sed -i 's/YYYYY/${ligand_name}.mol2/' $config
	sed -i 's/ZZZZZ/${receptor_name}_${ligand_name}/' $config
	PLANTS1.2_64bit --mode screen $config
	"""

}
