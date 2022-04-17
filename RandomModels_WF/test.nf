#!/usr/bin/env nextflow

params.foldx = "/home/nkhoza/generate_random_models/FoldX_ENCoM/testing/*.pdb"
params.encom = "path/to/*_Repair.pdb"
params.outfolder = "./FxE_output"

pdb = channel.fromPath(params.foldx)
process repair {

        input:
        file(pdb) from pdb

        output:
	file("*_Repair.pdb") into repaired        

        script:
        """
	ln -s /home/nkhoza/tools/FoldX/rotabase.txt 
	PATH=/home/nkhoza/tools/FoldX:$PATH
	foldx --command=RepairPDB --pdb=${pdb} 
        """
}

repaired.into {repaired_foldx; repaired_encom}

process foldx {

	input:
	file(repaired) from repaired_foldx

	script:
	"""
	ln -s /home/nkhoza/tools/FoldX/rotabase.txt
        PATH=/home/nkhoza/tools/FoldX:$PATH
    	foldx --command=Stability --pdb=${repaired}
	"""
}

process encom {

	input:
	file(repaired) from repaired_encom

	output:
	file("*.eigen") into DDGnDS
	
	script:
	"""
	PATH=/home/nkhoza/tools/ENCoM/bin:$PATH
	build_encom -i ${repaired} -cov ${repaired}.cov -o ${repaired}.eigen
	"""
}

process encom_processing {
	
	input:
	file(eigen) from DDGnDS
	
	script:
	"""
	sub_encom.py
	"""
}
