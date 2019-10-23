# ===============================
# Scripts to measure Radial Distribution function from the centre of mass
# of an arbitrary selection in a trajectory
# There will probably be issues if you try and use this on a trajectory with lots of atoms or frames
# vmd will probably run into memory issues if you try with trajectories > 1gb in size
# - Sam Wallace (wallace.samuel.j@gmail.com)
# ===============================

proc add_atoms_trajectory {molid atoms} {
# this is modified code from Axel Kohlmeyer's topotools mergemol process
# it will duplicate a trajectory by first copying a frame and all the bond
# information, then duplicating to match same number of frames as original
# trajectory. Once that is done, it runs another loop to set the coordinates 
# of all the atoms and the unit cell dimensions per frame. The extra atom(s) are 
# the last in the index.
package require topotools 1.1
    # compute total number of atoms and collect
    # offsets and number of atoms of each piece.
    set numatoms [molinfo $molid get numatoms]
	set numframes [molinfo $molid get numframes]
	set newtotal [expr {$numatoms + $atoms}]
    # create new molecule to hold data with same number of atoms and frames
    set mol -1
    if {[catch {mol new atoms $newtotal} mol]} {
        vmdcon -error "mergemols: could not create new molecule: $mol"
        return -1
    } else {
		set loopx 0
		while {$loopx < $numframes} {
			animate dup $mol
			set loopx [expr {$loopx + 1}]
			}
    }
    mol rename $mol [string range mergedmol-[join $molid -] 0 50]

    # copy data over piece by piece
    set bondlist {}
    set anglelist {}
    set dihedrallist {}
    set improperlist {}
	set oldsel [atomselect $molid all]
	set newsel [atomselect $mol "index 0 to [expr {$numatoms-$atoms}]"]

    # per atom props
    set cpylist {name type mass charge radius atomicnumber element x y z \
                         resname resid chain segname}
    $newsel set $cpylist [$oldsel get $cpylist]
		# assign structure data. we need to renumber indices
        set list [topo getbondlist both -molid $molid]
        foreach l $list {
            lassign $l a b t o
            lappend bondlist [list $a $b $t $o]
        }

        set list [topo getanglelist -molid $molid]
        foreach l $list {
            lassign $l t a b c 
            lappend anglelist [list $t $a $b $c]
        }

        set list [topo getdihedrallist -molid $molid]
        foreach l $list {
            lassign $l t a b c d
            lappend dihedrallist [list $t $a $b $c $d]
        }
        set list [topo getimproperlist -molid $molid]
        foreach l $list {
            lassign $l t a b c d
            lappend improperlist [list $t $a $b $c $d]
        }
        $oldsel delete
        $newsel delete
   
    # apply structure info
    topo setbondlist both -molid $mol $bondlist
    topo setanglelist -molid $mol $anglelist
    topo setdihedrallist -molid $mol $dihedrallist
    topo setimproperlist -molid $mol $improperlist

    # set box to be largest of the available boxes
    set amax 0.0
    set bmax 0.0
    set cmax 0.0
    foreach m $molid {
        lassign [molinfo $m get {a b c}] a b c
        if {$a > $amax} {set amax $a}
        if {$b > $bmax} {set bmax $b}
        if {$c > $cmax} {set cmax $c}
    }
    molinfo $mol set {a b c} [list $amax $bmax $cmax]
	# this part should iterate through all the frames and set the coordinates of each atom
	set loopy 0
	while {$loopy < $numframes} {
		puts $loopy
		animate goto $loopy
		set loopy [expr {$loopy + 1}]
		set oldsel [atomselect $molid all]
		set newsel [atomselect $mol "index 0 to [expr {$numatoms-$atoms}]"]
		# set each atoms coordinates
		set cpylist {x y z}
		$newsel set $cpylist [$oldsel get $cpylist]
		molinfo $mol set {a b c} [molinfo $molid get {a b c}]
		}
	$oldsel delete
	$newsel delete		
	mol addrep $mol
    return $mol
}


proc rdf_com {molid selectiontext} {
# this process will duplicate a and add an atom to a trajectory
# with coordinates that follow the center of mass of selection text
# example rdf_com top "segname LIG"

# collect some variables
set orig_num_atoms [molinfo $molid get numatoms]
set measure_molid [add_atoms_trajectory $molid 1]
set numframes [molinfo $molid get numframes]
# set the selections
set a [atomselect $orig_molid $selectiontext]
set temp_sel [atomselect $measure_molid "index = $orig_num_atoms"]
# the loop that sets the coordinates of the first dummy atom to the centre of mass
# of selection text for all frames
set loopy 0
	while {$loopy < $numframes} {
		puts $loopy
		animate goto $loopy
		set loopy [expr {$loopy + 1}]
		set b [measure inertia $a]
		set b_x [lindex [lindex $b 0] 0]
		set b_y [lindex [lindex $b 0] 1]
		set b_z [lindex [lindex $b 0] 2]
		# set each atoms coordinates
		set xyzlist [list [list COM DUM $b_x $b_y $b_z]]
		$temp_sel set {name type x y z} $xyzlist
	}
}
