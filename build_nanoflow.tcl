# This file is used to generate a nanoflow system like below
#
#	..->..@...||							                	 ||..*.....
#	..->......|| 								                 ||.....*..
#	..->.*....||=================================||........
#	..->...@...........*.............................@.....
#	..->.........@..............*...........@...........*..
#	.*->.......*.............@..||.......*.................
#	..->......@......*..........||.....*......@............
#	..->...*..||=================================||...@....
#	..->..@...||							                	 ||........
#	..->......||							                	 ||..*.....
#
#	-> stands for force;
#	@ and * are ions;
#	. is water molecule;
#	= and || are nanomaterials that can be graphene or boron nitride
#   TCL forces are needed to make molecules flow through the nanochannel (NAMD)
#   
#	Date 2020 Apr 30
#	Author: MM.Xue
#   Try using VMD

package require psfgen
package require nanotube
package require solvate
package require autoionize
package require pbctools
set tcl_precision 3

# proc to make the all in the middle of the tube
proc genmidwall {width r m seg name {t "armchair"} {nl 1} {bl 0.1418} {bt "C-C"}} {
	graphene -lx $width -ly $width -type $t -nlayers 1
	set sel [atomselect top all]
	$sel set segname $seg
	$sel moveby [vecinvert [measure center $sel]]
	$sel moveby "0 0 $m"
	set save [atomselect top "x*x+y*y < $r*$r-1  and y<=0"]
	$save writepsf $name.psf
	$save writepdb $name.pdb
	$sel delete
	$save delete
	mol delete top
	puts "============== mid wall build ================"
}

# proc to make walls on two ends of the nanotube
proc genwall {width r d seg name {t "armchair"} {nl 1} {bl 0.1418} {bt "C-C"} } {
	graphene -lx $width -ly $width -type $t -nlayers 1
    set sel [atomselect top all]
    $sel set segname $seg
    $sel moveby [vecinvert [measure center $sel]]
    $sel moveby "0 0 $d"
    set save [atomselect top "not (x*x+y*y <= $r*$r +1)"]
    $save writepsf $name.psf
    $save writepdb $name.pdb
    $sel delete
    $save delete
	mol delete top
	puts "============== wall build ================"
}

# proc to make nanotube
proc gentube {length c_n c_m seg name {t "armchair"} {nl 1} {bl 0.1418} {bt "C-C"}} {
	nanotube -l $length -n $c_n -m $c_m -cc $bl -ma $bt
	set sel [atomselect top all]
	$sel set segname $seg
	$sel moveby [vecinvert [measure center $sel]]
	$sel writepsf $name.psf
	$sel writepdb $name.pdb
	$sel delete
	mol delete top
	puts "============== tube build ================"
}


proc genmodel {tube_length tube_c_n tube_c_m tube_seg tube_name \
{mw_num 1} {mw_move 0} mw_seg mw_name \
ws_width ws_seg ws_name} {
	gentube $tube_length $tube_c_n $tube_c_m $tube_seg $tube_name
	mol new $tube_name.psf
	mol addfile $tube_name.pdb
	set sel [atomselect top all]
	set xyz [measure minmax $sel]
#	set dx [expr ([lindex [lindex $xyz 1] 0] - [lindex [lindex $xyz 0] 0])]
#	set dy [expr ([lindex [lindex $xyz 1] 1] - [lindex [lindex $xyz 0] 1])]
	set dz [expr ([lindex [lindex $xyz 1] 2] - [lindex [lindex $xyz 0] 2])]
#	set r [expr sqrt(($dx*$dx+$dy*$dy)/2.0)/2.0]
	set r [lindex [lindex $xyz 1] 0]
	set disleft [expr [lindex [lindex $xyz 0] 2] -1]
	set disright [expr [lindex [lindex $xyz 1] 2] +1]
	mol delete top
	$sel delete
	genwall $ws_width $r $disleft ${ws_seg}1 ${ws_name}1
	genwall $ws_width $r $disright ${ws_seg}2 ${ws_name}2
	resetpsf
	readpsf ${ws_name}1.psf
	coordpdb ${ws_name}1.pdb
	readpsf ${ws_name}2.psf
	coordpdb ${ws_name}2.pdb
	writepsf $ws_name.psf
	writepdb $ws_name.pdb

	set mw_width 4
	if {$mw_num == 1} {
		genmidwall $mw_width $r $mw_move $mw_seg $mw_name
	} elseif {$mw_num == 2} {
		genmidwall $mw_width $r $mw_move ${mw_seg}1 ${mw_name}1
		genmidwall $mw_width $r [expr -$mw_move] ${mw_seg}2 ${mw_name}2
		resetpsf	
		readpsf ${mw_name}1.psf
		coordpdb ${mw_name}1.pdb
		readpsf ${mw_name}2.psf
		coordpdb ${mw_name}2.pdb
		writepsf $mw_name.psf
		writepdb $mw_name.pdb
	}
	resetpsf 
	readpsf $tube_name.psf
	coordpdb $tube_name.pdb
	if {$mw_num !=0 } {
		readpsf $mw_name.psf
		coordpdb $mw_name.pdb
	}
	readpsf $ws_name.psf
	coordpdb $ws_name.pdb
	writepsf tube_walls.psf
	writepdb tube_walls.pdb
	puts "============== tube_walls build ================"
}

# proc to solvate the tube_walls model
proc gensol {sol_length ifion {cation_name "SOD"} {anion_name "CLA"} {ion_conc 0.6}} {
	mol new tube_walls.psf
	mol addfile tube_walls.pdb
	set all [atomselect top all]
	set box [measure minmax $all]
	set sol_width [expr [lindex [lindex $box 1] 0] - [lindex [lindex $box 0] 0]]
	set sol_width2 [expr [lindex [lindex $box 1] 1] - [lindex [lindex $box 0] 1]]
	set sol_left [lindex [lindex $box 0] 2]
	set sol_right [lindex [lindex $box 1] 2]
	set l1 "[expr -($sol_width/2)] [expr -($sol_width2/2)] [expr $sol_left-$sol_length]"
	set l2 "[expr ($sol_width/2)] [expr ($sol_width2/2)] [expr $sol_left-2]"
	set l3 "[expr -($sol_width/2)] [expr -($sol_width2/2)] [expr $sol_right+2]"
	set l4 "[expr ($sol_width/2)] [expr ($sol_width2/2)] [expr $sol_right+$sol_length]"
	solvate tube_walls.psf tube_walls.pdb -minmax "{$l1} {$l2}" -o tube_walls_sol_tmp  -s "WL"\
	-spsf ../../common/tip4p-kit/tip4p-20Abox.psf -spdb ../../common/tip4p-kit/tip4p-after1ns.pdb -stop ../../common/tip4p-kit/tip4p.top -ws 20 -ks "name OH2"
	solvate tube_walls_sol_tmp.psf tube_walls_sol_tmp.pdb -minmax "{$l3} {$l4}" -o tube_walls_sol -s "WR"\
	-spsf ../../common/tip4p-kit/tip4p-20Abox.psf -spdb ../../common/tip4p-kit/tip4p-after1ns.pdb -stop ../../common/tip4p-kit/tip4p.top -ws 20 -ks "name OH2"
	$all delete 
	mol delete top
	puts "============== solvated ================"
	if {$ifion == 1} {
		autoionize -psf tube_walls_sol.psf -pdb tube_walls_sol.pdb -sc $ion_conc -cation $cation_name -anion $anion_name -from 5 -between 5 -seg ION -o tube_walls_sol_ion
		mol new tube_walls_sol_ion.psf
		mol addfile tube_walls_sol_ion.pdb
		set carbon [atomselect top "name C"]
		set save [atomselect top all]
		$save set occupancy 0
		$carbon set occupancy 1
		$save writepsf system.psf
		$save writepdb system.pdb
		$save delete
		$carbon delete
		set watlist [atomselect top "name OH2"]
		set watlistid [$watlist get index]
		puts "[lindex $watlistid 0] and [expr [llength $watlistid]*3+[lindex $watlistid 0]]"
		$watlist delete
		mol delete top
		puts "============== ionized ================"
		puts "============== model sytem builded ================"
	} else {
		mol new tube_walls_sol.psf
		mol addfile tube_walls_sol.pdb
		set carbon [atomselect top "name C"]
		set save [atomselect top all]
		$save set occupancy 0
		$carbon set occupancy 1
		$save writepsf system.psf
		$save writepdb system.pdb
		$carbon delete
		$save delete
		set watlist [atomselect top "name OH2"]
		set watlistid [$watlist get index]
		puts "[lindex $watlistid 0] and [expr [llength $watlistid]*3+[lindex $watlistid 0]]"
		$watlist delete
		mol delete top
		puts "============== model system builded ================"
	}
}

#genmodel 8 8 8 "TB" "tube" 1 0 "MW" "mid_wall" 4 "WS" "walls" 
# var1: tube_length		var2: tube_c_m		var3: tube_c_n
# var4: midwall_num		var5: midwall_gap
# var6: sidwall_width
# var7: resorvoir_width	var8: ifion
genmodel 8 8 8 "TB" "tube" 0 0 "MW" "mid_wall" 4 "WS" "walls" 
gensol 75 1
