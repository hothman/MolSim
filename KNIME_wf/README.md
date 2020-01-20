# README

The KNIME workflow takes a multiple lines `sim` file and outputs a SDF file.

1- Read sim file
2- Strip salt molecule
3- Add Hydrogens 
4- Clean geometry
5- Generate conformers
6- Output to `SDF` file

Requires RDKit and CDK (Chemistry Development Kit) nodes. 
