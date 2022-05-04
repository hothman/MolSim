#!/usr/bin/env nextflow

params.foldx = "/media/houcem/theDrum/BILIM/nhlamulo/MolSim/RandomModels_WF/test/*.pdb"
params.rotabase = "/media/houcem/theDrum/modules/foldx/rotabase.txt"
params.output= "entropy_ddG_all.csv"

pdb = channel.fromPath(params.foldx)
process repair {

        input:
        file(pdb) from pdb

        output:
	file("*_Repair.pdb") into repaired        

        script:
        """
        ln -s ${params.rotabase}
	foldx --command=RepairPDB --pdb=${pdb} 
        """
}

repaired.into {repaired_foldx; repaired_encom}

process foldx {

	input:
	file(repaired) from repaired_foldx
	
	output:
	file("*_ST.fxout") into dG

	script:
	"""
	ln -s ${params.rotabase}
    	foldx --command=Stability --pdb=${repaired}
	"""
}

process encom {

	input:
	file(repaired) from repaired_encom

	output:
	file("*.cov") into cov
	val(name) into the_basename
	
	script:
		name = repaired.baseName.replaceFirst(".pdb","")

	"""
	build_encom -i ${repaired} -cov ${name}.cov -o ${name}.eigen
	"""
}

process processOutputs {
	echo true
	input:
	file(cov_file) from cov
	file(dg_file) from dG
	val(name) from the_basename
	
	output:
	file("*.csv") into dS
	
	script:
	"""
	get_entropy_dG.py --cov ${cov_file} --foldx ${dg_file} --output ${name}.csv

	"""
}

process outputCsv {
	echo true
	publishDir "./", mode:'copy'
	input: 
		file(all_data) from dS.collect()
	output: 
		file(params.output)
	 
	"""
	echo "tag_encom,S_vib,tag_foldx,dG" > ${params.output}
	cat $all_data >> ${params.output}
	"""
}
