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
	for pdb in ${pdb}
	do
		foldx --command=RepairPDB --pdb=\$pdb
	done 
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
	for repaired in ${repaired}
	do
    	    foldx --command=Stability --pdb=\$repaired
    	done
	"""
}

process encom {

	input:
	file(repaired) from repaired_encom

	output:
	
	script:
	"""
	PATH=/home/nkhoza/tools/ENCoM/bin:$PATH
	for repaired in ${repaired}
	do
		build_encom -i \$repaired -cov \$repaired.cov -o \$repaired.eigen
	done
	"""
}