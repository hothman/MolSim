# echo $2 $3 | gmx pdb2gmx -f $1  -o peptide.gro

gmx pdb2gmx -f 17aa.pdb -o peptide.gro -ss yes  << EOF
3
1
y
EOF

cat <<EOF  >merge_groups.ndx
q
EOF

cat << EOF > ions.mdp
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator	= steep		; Algorithm (steep = steepest descent minimization)
emtol		= 1000.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01      ; Energy step size
nsteps		= 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
coulombtype	    = PME		; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		; Short-range electrostatic cut-off
rvdw		    = 1.0		; Short-range Van der Waals cut-off
pbc		        = xyz 		; Periodic Boundary Conditions (yes/no)

EOF

gmx editconf -f peptide.gro -o boxed.gro -c -d 1.2 -bt octahedron
gmx solvate -cs -cp boxed.gro -o solvated.gro -p topol.top
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2
gmx genion -s ions.tpr -o solvated.gro -p topol.top -pname NA -nname CL -neutral << EOF
12
EOF

cat << EOF >em1.mdp
------em.mdp------
integrator = steep
nsteps = 1000
nstlist = 10
rlist = 1.0
coulombtype = pme
rcoulomb = 1.0
vdw-type = cut-off
rvdw = 1.0
nstenergy = 10
define = -DPOSRES
------------------
EOF

cat << EOF >em2.mdp
------em.mdp------
integrator = cg
nsteps = 200
nstlist = 10
rlist = 1.0
coulombtype = pme
rcoulomb = 1.0
vdw-type = cut-off
rvdw = 1.0
nstenergy = 10
------------------
EOF

gmx grompp -f em1.mdp -c solvated.gro -p topol.top -o em1.tpr
gmx mdrun -v -deffnm em1

TEMP=298
num_cpu=8
GRO_START=em2.gro

cat << EOF >nvt_equi_temperature.mdp

title		= OPLS Lysozyme NVT equilibration 
define		= -DPOSRES	; position restrain the protein
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= 10000		; 2 ps as in Elio et al (2012)
dt		    = 0.002		; 2 fs
; Output control
nstxout		= 1000		; save coordinates every 1.0 ps
nstvout		= 1000		; save velocities every 1.0 ps
nstenergy	= 1000		; save energies every 1.0 ps
nstlog		= 1000		; update log file every 1.0 ps
; Bond parameters
continuation	        = no		; first dynamics run
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10		; 20 fs, largely irrelevant with Verlet
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME	; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		; cubic interpolation
fourierspacing	= 0.12	; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1           ; time constant, in ps
ref_t		= $TEMP 	 $TEMP           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl		= no 		; no pressure coupling in NVT
; Periodic boundary conditions
pbc		= xyz		    ; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= yes		; assign velocities from Maxwell distribution
gen_temp	= 310		; temperature for Maxwell distribution
gen_seed	= -1		; generate a random seed

EOF

gmx grompp -f nvt_equi_temperature.mdp -p topol.top -c em1.gro -o equi1.tpr -maxwarn 2

gmx mdrun -pin on  -nt $num_cpu -v -deffnm equi1



GRO_START=equi1.gro
CHECKPOINT=equi1.cpt
TOPOL=topol.top
STEPNUM=150000000
num_cpu=8
TEMP=298
GROUPS_NDX=groups.ndx

gmx make_ndx -f $GRO_START -o $GROUPS_NDX <merge_groups.ndx

cat << EOF >production_run.mdp 

title       = Protein-ligand complex MD simulation 
; Run parameters
integrator  = md        ; leap-frog integrator
nsteps      = $STEPNUM    ; 2 * 1500000 = 1000 ps (3 ns)
dt          = 0.002     ; 2 fs
; Output control
;nstxout             = 0         ; suppress .trr output 
;nstvout             = 0         ; suppress .trr output
nstenergy           = 1000      ; save energies every 10.0 ps
nstlog              = 1000      ; update log file every 10.0 ps
nstxout-compressed  = 2500      ; write .xtc trajectory every 5.0 ps
compressed-x-grps   = System
energygrps          = Protein 
; Bond parameters
continuation    = yes           ; first dynamics run
constraint_algorithm = lincs    ; holonomic constraints 
constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
lincs_iter      = 1             ; accuracy of LINCS
lincs_order     = 4             ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type         = grid      ; search neighboring grid cells
nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb        = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw            = 1.0       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order       = 4         ; cubic interpolation
fourierspacing  = 0.12      ; grid spacing for FFT
; Temperature coupling
tcoupl      = V-rescale                     ; modified Berendsen thermostat
tc-grps     = Protein non-Protein    ; two coupling groups - more accurate
tau_t       = 0.1    0.1                    ; time constant, in ps
ref_t       = $TEMP   $TEMP                     ; reference temperature, one for each group, in K
; Pressure coupling 
pcoupl      = Parrinello-Rahman             ; pressure coupling is on for NPT
pcoupltype  = isotropic                     ; uniform scaling of box vectors
tau_p       = 2.0                           ; time constant, in ps
ref_p       = 1.0                           ; reference pressure, in bar
compressibility = 4.5e-5                    ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc         = xyz       ; 3-D PBC
; Dispersion correction
DispCorr    = EnerPres  ; account for cut-off vdW scheme
; Velocity generation
gen_vel     = no        ; assign velocities from Maxwell distribution


EOF

cat $GRO_START > temp_coor.gro
cat $CHECKPOINT > temp_checkpoint.cpt

i="1"
while [ $i -le 1 ]

do 


gmx grompp -f production_run.mdp -c temp_coor.gro -t  temp_checkpoint.cpt -p $TOPOL -n $GROUPS_NDX -o md_$i.tpr 


echo "##########################"
echo "      production $i       "
echo "##########################"


gmx mdrun -nt $num_cpu -pin on  -deffnm md_$i -o md_$i.trr

echo 1 1  | gmx trjconv -center -f md_$i.xtc -n groups.ndx -o  md_nowater_$i.xtc -s md_$i.tpr -pbc mol -ur compact 

rm md_$i.xtc
cat md_$i.gro >temp_coor.gro
cat md_$i.cpt >temp_checkpoint.cpt 

i=`expr $i + 1`

done

echo 1 1  | gmx trjconv -center -f $GRO_START  -o  complex_no_water.gro -s md_1.tpr -pbc mol -ur compact

