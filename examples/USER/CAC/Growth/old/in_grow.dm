	units       metal

	dimension    3
	boundary     p p f
	comm_style    cac
	newton off

#Manually define processor subdomains 
	processors 16 12 1

	#read_data     Pb.txt
	read_restart   restart.50000

	change_box all z final -0.5 950

	pair_style cac/pb/dsf 12.0 0.24 12.5 one
	pair_coeff 1 2 4880000 0.1730 211.0
	pair_coeff 1 3 8434000 0.1765 241.0
	pair_coeff 2 2 5200.00 0.3840 127.0
	pair_coeff 3 3 4475.00 0.3825 130.0
	pair_coeff 2 3 5060.00 0.3860 129.0

	interface_quadrature off

	timestep     0.01
	atom_modify first atom

	variable h equal 606	#position of numerical interface (Atoms/Elements)
	variable rb equal 940	#position of reflective boundary

	region		1 block EDGE EDGE EDGE EDGE EDGE 9.5 units box	#Bottom layer of substrate
	region		2 block EDGE EDGE EDGE EDGE 9.5 $(v_h-1) units box	#CG region
    region		3 block EDGE EDGE EDGE EDGE $(v_h+5) $(v_h+15) units box	# Heat region
	region		4 block EDGE EDGE EDGE EDGE $(v_h+25) $(v_h+50) units box	#Region to check temperature
	region      5 block EDGE EDGE EDGE EDGE $(v_h-1) EDGE units box	#Atom region
	region      6 block EDGE EDGE EDGE EDGE 9.5 375.5 units box #CG region: Element size 16
	region      7 block EDGE EDGE EDGE EDGE 375.5 558.5 units box #CG region: Element size 8
	region      8 block EDGE EDGE EDGE EDGE 558.5 605 units box	#CG region: Element size 4
	
	# Region 9: Heat region when only run NVE on the atoms at substrate top surface
	region		9 block EDGE EDGE EDGE EDGE $(v_h+15) $(v_h+22) units box	
	# Region 10: The atoms at substrate top surface
	region		10 block EDGE EDGE EDGE EDGE $(v_h+15) EDGE units box
	
	region      20 block EDGE EDGE EDGE EDGE $(v_rb-50) $(v_rb) units box  #region for atom creation
	region		21 block EDGE EDGE EDGE EDGE $(v_rb-52) EDGE units box 	 
	delete_atoms region 21	#Delete the atoms in the region

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
	group       CG4 region 9
	group       CG4 region 10

	compute cgtemp CG cac/nodal/temp
	compute cgtemp4 CG4 cac/nodal/temp
	compute cgtemp8 CG8 cac/nodal/temp
	compute cgtemp16 CG16 cac/nodal/temp
	compute subtemp1 heatD1 cac/nodal/temp
	compute subtemp2 heatD2 cac/nodal/temp
	thermo_style custom step c_cgtemp c_cgtemp4 c_cgtemp8 c_cgtemp16 c_subtemp1 c_subtemp2
	thermo_modify lost ignore flush yes
	thermo	     200

	restart 100000 restart.*

#Define processor subdomains using rcb
#required to weight the elements and load balance non-uniform #meshes
	#compute Eweight all cac/quad/count
	#variable Eweights atom c_Eweight
	#fix comm all balance 1 1.00 rcb weight var Eweights

#Define processor subdomains usingshift
	#fix comm all balance 1000 1.1 shift z 5 1.05 weight var Eweights

	#dump	     1 all cac/nodal/positions 1000000 CACmesh.txt
	#dump         2 atom cac/xyz 20000 dump.xyz
	dump    	 3 all custom 2000 dump_all.xyz id type x y z

	fix	   reflect all cac/oneway 1 21 -z
	fix    3 bottom cac/setforce 0.0 0.0 0.0
	fix    4 bottom cac/setvelocity 0.0 0.0 0.0

	fix Momentum1 CG8 cac/momentum 10 linear 1 1 0 angular 	# Avoid rotation
	fix Momentum2 CG16 cac/momentum 10 linear 1 1 0 angular
	run 0
	
#########################################################################################
	#LOOP
	variable a loop 1000
	label loopa

	create_atoms 1 random 400 $(v_a*10+13) 20
	create_atoms 2 random 400 $(v_a*20+13) 20
	
	set type 1 cac/charge 0.8
	set type 2 cac/charge -0.8

	group       free subtract all bottom
	group       atom region 5
	group       heat_top region 9
	group       free_top region 10
	group 		rainD dynamic all region 21 every 1

	fix    5 rainD cac/setvelocity 0.0 0.0 -1.0

#########################################################################################
	#unfix comm
	#fix comm all balance 200 1.1 shift z 5 1.05 weight var Eweights

	comm_modify group all
	neigh_modify every 1 delay 5 check yes include all
	fix NVE free cac/nve/limit 0.5

	fix TCA heat1 cac/viscous/temp 0.020 910 20

	fix TC1 CG4 cac/viscous/temp 0.03 910.0 5
	fix TC2 CG8 cac/viscous/temp 0.03 55.0 5
	fix TC3 CG16 cac/viscous/temp 0.03 3.0 5

	fix TF1 CG8 cac/tempforce 18.8524 21.7707 18.8530 65.53598 900.0 8.0
	fix TF2 CG16 cac/tempforce 18.8524 21.7707 18.8530 65.53598 900.0 8.0
	run 5000

	unfix NVE
	unfix TCA
	unfix TC1
	unfix TC2
	unfix TC3
	unfix TF1
	unfix TF2
#########################################################################################
	#unfix comm
	#fix comm all balance 2000 1.1 shift z 5 1.05 weight group 3 bottom 1 CG 1 atom 3

	comm_modify cutoff 12.5 group atom
	neigh_modify every 1 delay 5 check yes include atom
	fix NVE free_top cac/nve/limit 0.5
	fix TCA heat_top cac/viscous/temp 0.01 910 20
	run 45000
	unfix NVE
	unfix TCA

	unfix 5
#########################################################################################
	delete_atoms region 21	#Delete the atoms in the region
	next a
	jump in_grow.dm loopa