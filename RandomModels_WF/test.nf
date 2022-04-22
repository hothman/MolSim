#!/usr/bin/env nextflow

params.foldx = "/path_to/*.pdb"
params.rotabase = "/path/to/rotabase/rotabase.txt"
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

process encom_processing {
	echo true
	input:
	file(cov_file) from cov
	val(name) from the_basename
	
	output:
<<<<<<< HEAD
	file("*.csv") into dS
	
	script:
	"""
	get_entropy.py --cov ${cov_file} --output ${name}.csv
=======
	file("${cov}_output.cov") into dS
	
	script:
	"""
	get_entropy.py --cov ${cov} --output ./${cov}_output.csv
>>>>>>> 8294394c7f3c9616becade037827b8eea33d7529
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
	echo "tag,S_vib" > ${params.output}
	cat $all_data >> ${params.output}
	"""
}
