#set topname [glob *top]
#foreach i $topname {
#    if { [catch {cg_bonds -gmx /usr/local/gromacs-5.0.7/bin/gmxdump -top $i} errmsg] } {
#    } else {
#        break
#    }
#}
set tprname [glob *tpr]
foreach i $tprname {
    if { [catch {cg_bonds -gmx /home/xuemm/software/gromacs-5.0.7/bin/gmxdump -tpr $i -cutoff 8.0 } errmsg] } {
    } else {
        puts "cg_bonds -gmx /usr/local/gromacs-5.0.7/bin/gmxdump -tpr $i"
        break
    }
}

package require topotools
package require pbctools
#new material for SMA copolymer
#material add SMA copy EdgyShiny
#material change {ambient specular diffuse shininess mirror opacity outline outlinewidth} SMA {0.05 0.80 0.70 0.70 0.00 1.00 0.70 0.80}
#Display settings
#color Display Background white
pbc box -on
#display projection orthographic
#display rendermode GLSL
#display shadows on
#display ambientocclusion on
#display aoambient 0.90
#display cuedensity 0.30
#display height 2.0

mol delrep 0 top
set options {lipid po4 nc3 sma carbonxyl}
array set colors {lipid "ColorID 8" po4 "ColorID 3" nc3 "ColorID 0" sma "ColorID 7" carbonxyl "ColorID 4"}
array set styles {lipid "Licorice 0.5 15 15" po4 "VDW 1.2 15" nc3 "VDW 1.2 15" sma "Licorice 0.5 15 15" carbonxyl "VDW 1.2 15"}
array set materi {lipid "AOChalky" po4 "AOChalky" nc3 "AOChalky" sma "SMA" carbonxyl "SMA"}
array set select {lipid "not resname SM1 SM2 SM3 SM4 SM5 PW ION and same resid as name PO4" po4 "name PO4" nc3 "name NC3" sma "resname SM1 SM2 SM3 SM4 SM5" carbonxyl "name CB1 CB2"}

foreach idoption $options {
    mol selection $select($idoption)
    mol rep $styles($idoption)
    mol material $materi($idoption)
    mol color $colors($idoption)
    mol addrep top
}
    
