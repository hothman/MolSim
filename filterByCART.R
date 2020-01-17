# Author Houcemeddine Othman
# Wits University 2019 

# Implementation of the CART model to filter ligands with 
# high potency to cross the BBB, Check Suenderhauf et al, 2012 (PMC6269008)

library("Rcpi")
library('tidyverse')
Input = 'AllCompounds_ChemblPAINS.smi'

smiles <- read.table(Input)
smi <- readMolFromSmi(Input, type = "mol")
BCUT = extractDrugBCUT(smi)
 #= dat$BCUTw.1l
TPSA = extractDrugTPSA(smi)
aLogP = extractDrugALOGP(smi)


my_descriptors <- data.frame(smiles$V1, BCUT$BCUTw.1l, TPSA$TopoPSA, aLogP$ALogP )
colnames(my_descriptors ) <- c( "smiles", "BCUT", "TPSA", "aLogP" )

filtered <- filter(my_descriptors, (my_descriptors$aLogP <= -0.4263 & my_descriptors$BCUT <= 11.9)  |  
         (my_descriptors$aLogP > -0.4263 & my_descriptors$TPSA <= 150.54553))

print.data.table(filtered)
write.table(filtered$smiles, 'AllCompounds_ChemblPAINSBBB.smi', quote=FALSE, row.names=FALSE, col.names=FALSE )
