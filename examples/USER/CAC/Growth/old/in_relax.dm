units       metal

dimension    3
boundary     p p f
comm_style    cac
atom_style     cac/charge 8 2
newton off

#Manually define processor subdomains 
processors 16 12 1

read_data     Pb.txt
#read_restart   restart.7250000

pair_style cac/pb/dsf 12.0 0.24 12.5 one
pair_coeff 1 2 4880000 0.1730 211.0
pair_coeff 1 3 8434000 0.1765 241.0
pair_coeff 2 2 5200.00 0.3840 127.0
pair_coeff 3 3 4475.00 0.3825 130.0
pair_coeff 2 3 5060.00 0.3860 129.0

interface_quadrature off

timestep     0.01
atom_modify first atom

	region		1 block EDGE EDGE EDGE EDGE EDGE 9.5 units box	#Bottom layer of substrate
	region		2 block EDGE EDGE EDGE EDGE 9.5 605 units box	#CG region
    region		3 block EDGE EDGE EDGE EDGE 605 616 units box	# Heat region
	region		4 block EDGE EDGE EDGE EDGE 625 635 units box	#Region to check temperature
	region      5 block EDGE EDGE EDGE EDGE 606 EDGE units box	#Atom region
	region      6 block EDGE EDGE EDGE EDGE 9.5 375.5 units box #CG region: Element size 16
	region      7 block EDGE EDGE EDGE EDGE 375.5 558.5 units box #CG region: Element size 8
	region      8 block EDGE EDGE EDGE EDGE 558.5 605 units box	#CG region: Element size 4
	
	region      11 block EDGE EDGE EDGE EDGE 650 EDGE units box #Reflect boundary
	
	group       bottom region 1
	group       free subtract all bottom
	group       CG region 2
	group       heatD1 dynamic all region 3 every 100
	group       heat1 region 3
	group       heatD2 dynamic all region 4 every 100
	group       heat2 region 4
	group       atom region 5	
	group       CG16 region 6
	group       CG8 region 7
	group       CG4 region 8


compute cgtemp CG cac/nodal/temp
compute cgtemp4 CG4 cac/nodal/temp
compute cgtemp8 CG8 cac/nodal/temp
compute cgtemp16 CG16 cac/nodal/temp
compute subtemp1 heatD1 cac/nodal/temp
compute subtemp2 heatD2 cac/nodal/temp
thermo_style custom step c_cgtemp c_cgtemp4 c_cgtemp8 c_cgtemp16 c_subtemp1 c_subtemp2
thermo_modify lost ignore flush yes
thermo	     100

#Define processor subdomains using rcb
#required to weight the elements and load balance non-uniform #meshes
#compute Eweight all cac/quad/count
#variable Eweights atom c_Eweight
#fix comm all balance 1 1.00 rcb weight var Eweights

#Define processor subdomains usingshift
#fix comm all balance 1000 1.1 shift z 5 1.05 weight var Eweights

restart 10000 restart.*

#dump	     1 all cac/nodal/positions 1000000 CACmesh.txt #Output dump file for Tecplot (Need post-processing)
#dump        2 atom cac/xyz 50000 dump_all.xyz	#Output dump file for OVITO, atoms in elements are output as well, file will be huge
dump    	 3 all custom 2000 dump.xyz id type x y z #Output dump file for OVITO


fix	   reflect all cac/oneway 1 11 -z
fix    3 bottom cac/setforce 0.0 0.0 0.0
fix    4 bottom cac/setvelocity 0.0 0.0 0.0

run 0

fix Momentum1 CG8 cac/momentum 5 linear 1 1 0 angular	# Avoid rotation
fix Momentum2 CG16 cac/momentum 5 linear 1 1 0 angular

#INITIAL TEMP
fix NVE free cac/nve
#fix TC_AT heat1 cac/temp/rescale 20 900 900 5 0.95

#Temperature control of CG elements	
fix TC1 CG4 cac/viscous/temp 0.03 950.0 5 
fix TC2 CG8 cac/viscous/temp 0.03 100.0 5
fix TC3 CG16 cac/viscous/temp 0.03 5.0 5

#Forces for thermal expansion, the last factor requires careful testing, to match thermal expansion rate
fix TF1 CG8 cac/tempforce 18.8524 21.7707 18.8530 65.53598 900.0 8.0 
fix TF2 CG16 cac/tempforce 18.8524 21.7707 18.8530 65.53598 900.0 8.0

run 30000
fix TC_AT heat1 cac/temp/rescale 20 900 900 5 0.1	#Heat to growth temperature of 900K
run 5000
unfix TC_AT

fix TC heat1 cac/viscous/temp 0.05 900 10 #Relax the substrate at the growth temperature of 900K
run 15000

