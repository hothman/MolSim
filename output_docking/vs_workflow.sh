#!/bin/bash 

# smiles file 
SMILES=./data/substances.smi 

while read -r line 
do
	echo $line >lig.smi 
	# get Dir name as ZINC identifier
	dirname=$(awk {'print $2'} lig.smi)
	echo "######### Docking ligand $dirname ###############"
	obabel -i smi lig.smi -o pdbqt -O lig.pdbqt
	mkdir $dirname 
	mv lig.smi $dirname
	mv lig.pdbqt $dirname
	
done<$SMILES 
