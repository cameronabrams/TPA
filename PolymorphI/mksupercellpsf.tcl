package require pbctools
mol new my_ub7_unitcell.psf
mol addfile my_ub7_unitcell.pdb
set base [atomselect top all]
set resids [lsort -unique [$base get resid]]
set segids [lsort -unique [$base get segname]]
set chains [lsort -unique [$base get chain]]
puts "$resids $segids $chains"
set nr_c [llength $resids]
set na_c [$base num]
set ns_c [llength $segids]
set nc_c [llength $chains]

set d2r [expr $M_PI/180.0]
set p [pbc get]
set a [lindex [lindex $p 0] 0]
set b [lindex [lindex $p 0] 1]
set c [lindex [lindex $p 0] 2]
set alpha [expr [lindex [lindex $p 0] 3]*$d2r]
set beta [expr [lindex [lindex $p 0] 4]*$d2r]
set gamma [expr [lindex [lindex $p 0] 5]*$d2r]

set ucellbasisvector1 [list $a 0 0]
set ucellbasisvector2 [list [expr $b*cos($gamma)] [expr $b*sin($gamma)] 0]
set ucellbasisvector3 [list [expr $c*cos($beta)] [expr $c*(cos($alpha)-cos($beta)*cos($gamma))/sin($gamma)] [expr $c*sqrt(1-pow(cos($beta),2)-pow((cos($alpha)-cos($beta)*cos($gamma))/sin($gamma),2))]]

set ucellV [expr $a*$b*$c*sqrt(1+2*cos($alpha)*cos($beta)*cos($gamma)-pow(cos($alpha),2)-pow(cos($beta),2)-pow(cos($gamma),2))]

set x0 [$base get {x y z}]
set r0 [$base get resid]

set nx 8
set ny 6
set nz 12

set nc [expr $nx * $ny * $nz]

set sc_na [expr $nc*$na_c]
set sc_nr [expr $nc*$nr_c]

set sc_a [expr $nx * $a]
set sc_b [expr $ny * $b]
set sc_c [expr $nz * $c]

puts "Supercell $sc_a x $sc_b x $sc_c $nc cells with $sc_na atoms in $sc_nr residues"

set fp [open "cell.inp" "w"]
puts $fp "cellbasisvector1 [expr $nx*[lindex $ucellbasisvector1 0]] [expr $ny*[lindex $ucellbasisvector1 1]] [expr $nz*[lindex $ucellbasisvector1 2]]"
puts $fp "cellbasisvector2 [expr $nx*[lindex $ucellbasisvector2 0]] [expr $ny*[lindex $ucellbasisvector2 1]] [expr $nz*[lindex $ucellbasisvector2 2]]"
puts $fp "cellbasisvector3 [expr $nx*[lindex $ucellbasisvector3 0]] [expr $ny*[lindex $ucellbasisvector3 1]] [expr $nz*[lindex $ucellbasisvector3 2]]"
puts $fp "cellorigin 0 0 0"
close $fp

set sc_cryst [list $sc_a $sc_b $sc_c [expr $alpha/$d2r] [expr $beta/$d2r] [expr $gamma/$d2r]]

set ci 0
for {set i 0} {$i < $nx} {incr i} {
    set imvb [vecscale $ucellbasisvector1 $i]
    for {set j 0} {$j < $ny} {incr j} {
        set jmvb [vecscale $ucellbasisvector2 $j]
        for {set k 0} {$k < $nz} {incr k} {
            set kmvb [vecscale $ucellbasisvector3 $k]
            set mvb [vecadd $imvb [vecadd $jmvb $kmvb]]
            $base set {x y z} $x0
            $base set resid $r0
            $base moveby $mvb
            set this_resid []
            for {set r 0} {$r < [llength $r0]} {incr r} {
                lappend this_resid [expr [lindex $r0 $r]+2*$ci]
            }
            $base set resid $this_resid
            $base writepdb "${ci}.pdb"
            incr ci
        }
    }
}

package require psfgen
topology /home/cfa/charmm/toppar/top_all36_prot.rtf
topology /home/cfa/charmm/toppar/top_all36_carb.rtf
topology /home/cfa/charmm/toppar/top_all36_lipid.rtf
topology /home/cfa/charmm/toppar/top_all36_na.rtf
topology /home/cfa/charmm/toppar/top_all36_cgenff.rtf
topology ../charmm/ub7.str

segment A {
    set ci 0
    for {set i 0} {$i < $nx} {incr i} {
        for {set j 0} {$j < $ny} {incr j} {
            for {set k 0} {$k < $nz} {incr k} {
                pdb ${ci}.pdb
                incr ci
            }
        }
    }
}
set ci 0
for {set i 0} {$i < $nx} {incr i} {
    for {set j 0} {$j < $ny} {incr j} {
        for {set k 0} {$k < $nz} {incr k} {
            coordpdb ${ci}.pdb A
            incr ci
#            exec rm ${ci}.pdb
        }
    }
}
writepsf my_ub7_supercell.psf
writepdb my_ub7_supercell_tmp.pdb

mol new my_ub7_supercell_tmp.pdb
pbc set $sc_cryst
[atomselect top all] writepdb my_ub7_supercell.pdb
exit
