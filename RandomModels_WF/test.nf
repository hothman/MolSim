#!/usr/bin/env nextflow

params.foldx = "/home/nkhoza/generate_random_models/FoldX_ENCoM/testing/*.pdb"

pdb = channel.fromPath(params.foldx)
process repair {

        input:
        file(pdb) from pdb

        output:
	file("*_Repair.pdb") into repaired        

        script:
        """
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
    	foldx --command=Stability --pdb=${repaired}
	"""
}

process encom {

	input:
	file(repaired) from repaired_encom

	output:
	file("*.cov") into cov
	
	script:
	"""
	build_encom -i ${repaired} -cov ${repaired}.cov -o ${repaired}.eigen
	"""
}

process encom_processing {
	
	input:
	file(cov) from cov
	
	output:
	file(*_output.cov) into	dS
	
	script:
	"""
	get_entropy.py --cov ${cov} --output ./{cov)_output.csv
	"""
}

entopy = dS.collectFile(name: "./entropy_all.csv", newLine: true, skip: 1, keepHeader: true)
folding_energy = dG.collectFile(name: "./folding_energy_all.csv", newLine: true, skip: 1, keepHeader: true)
