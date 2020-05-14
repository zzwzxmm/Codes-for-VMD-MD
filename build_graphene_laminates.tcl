# This script is used to build a model system of graphene laminates merged in ion solvation
# Note the graphene sheets here are of charges
# Modify before using 
# Author: MINMIN XUE
# Cite this article if using this script:
#               Heterogeneous graphene oxide membrane for rectified ion transport
#               Wenwen Fei, Minmin Xue, Hu Qiu, Wanlin Guo,2019,Nanoscale


package require psfgen
package require nanotube
package require solvate
package require autoionize
package require pbctools
set tcl_precision 3

proc sheet { seg_name charge_value lx ly type name_single dir} {
    graphene -lx [expr $lx/10] -ly [expr $ly/10] -type $type -nlayers 1
   	set sel [atomselect top all]
   	$sel set charge $charge_value
	$sel set segname $seg_name
	$sel writepsf $dir/$name_single.psf
	$sel writepdb $dir/$name_single.pdb
	mol delete top
	$sel delete
}

proc duplicate_sheet { name_single lx ly gap_x gap_y gap_z name layer_number sheet_number gap_move dir} {
	mol new $dir/$name_single.psf
	mol addfile $dir/$name_single.pdb
	set sel [atomselect top all]
	set seg_name [lindex [$sel get segname] 0]

	for {set j 0} {$j < $layer_number} {incr j} {
		for {set i 0} {$i < $sheet_number} {incr i} {
			set a [expr ($gap_x+$lx)*$i+$gap_move*($j%2)]
			set b 0
			set c [expr $gap_z*$j]
			set vec "$a $b $c"
			set vec_inv [vecinvert $vec]
			$sel moveby $vec
			$sel set segname $seg_name$i$j
			$sel writepsf $dir/$name$i$j.psf
			$sel writepdb $dir/$name$i$j.pdb
			$sel moveby $vec_inv
		}
	}
	$sel delete

	resetpsf 
	for {set j 0} {$j < $layer_number} {incr j} {
		for {set i 0} {$i < $sheet_number} {incr i} {
			readpsf $dir/$name$i$j.psf
			coordpdb $dir/$name$i$j.pdb
		}
	}
	writepsf $dir/$name.psf
	writepdb $dir/$name.pdb
	mol delete all
}

set lx_g 20; 			set lx_p 20
set ly_g 20; 			set ly_p 20
set type "armchair"; 		set type_2 "zigzag"
set pcharge 0.01; 		set ncharge -0.01
set seg_name_g "G";		set seg_name_p "P"
set name_single_g "gosheet"
set name_single_p "peisheet"
set gap_x_g 10 ;		set gap_x_p 10
set gap_y_g 6 ;			set gap_y_p 6
set gap_z_g 10 ;		set gap_z_p 9
set name_g "go";		set name_p "pei"
set layer_number_g 4;		set layer_number_p 4
set sheet_number_g 3;		set sheet_number_p 3
set gap_move_g 10;		set gap_move_p 10
set name "peigo";		set dir "pdb_psf_file"

sheet $seg_name_g $pcharge $lx_g $ly_g $type $name_single_g $dir
sheet $seg_name_p $ncharge $lx_p $ly_p $type $name_single_p $dir
duplicate_sheet $name_single_g $lx_g $ly_g $gap_x_g $gap_y_g $gap_z_g $name_g $layer_number_g $sheet_number_g $gap_move_g $dir
duplicate_sheet $name_single_p $lx_p $ly_p $gap_x_p $gap_y_p $gap_z_p $name_p $layer_number_p $sheet_number_p $gap_move_p $dir

set gap_pg "10"
mol new $dir/$name_p.psf
mol addfile $dir/$name_p.pdb
set sel [atomselect top all]
set gap_pg_z [expr $layer_number_g*$gap_z_g+$gap_pg]
$sel moveby "0 0 $gap_pg_z"
$sel writepsf $dir/$name_p.psf
$sel writepdb $dir/$name_p.pdb
mol delete top

resetpsf
readpsf $dir/$name_g.psf
coordpdb $dir/$name_g.pdb
readpsf $dir/$name_p.psf
coordpdb $dir/$name_p.pdb
writepsf $dir/$name.psf
writepdb $dir/$name.pdb
mol new $dir/$name.psf 
mol addfile $dir/$name.pdb
set pbc_a [expr ($lx_g+$gap_x_g)*$sheet_number_g]
set pbc_b [expr $ly_g+$gap_y_g]
set pbc_c [expr $gap_z_g*$layer_number_g+$gap_z_p*$layer_number_p+$gap_pg*4]
#puts "$a $b $c"
pbc set "$pbc_a $pbc_b $pbc_c"
set sel [atomselect top all]
$sel moveby "0 [expr $gap_y_g/2] [expr $gap_pg*3/2]"
$sel writepsf $dir/$name.psf
$sel writepdb $dir/$name.pdb 
mol delete top
set origin "0 0 0"
set maxv "$pbc_a $pbc_b $pbc_c"
set minmax "{$origin} {$maxv}"
solvate $dir/$name.psf $dir/$name.pdb -minmax $minmax  -o $name-sol -b 2.0
autoionize -psf $name-sol.psf -pdb $name-sol.pdb -sc 0.3 -cation POT -anion CLA -from 5 -between 5 -seg ION -o $name-ionize
