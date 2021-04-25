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

set p [pbc get]
set a [lindex [lindex $p 0] 0]
set b [lindex [lindex $p 0] 1]
set c [lindex [lindex $p 0] 2]

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
puts $fp "cellbasisvector1 $sc_a 0 0"
puts $fp "cellbasisvector2 0 $sc_b 0"
puts $fp "cellbasisvector3 0 0 $sc_c"
puts $fp "cellorigin 0 0 0"
close $fp

set sc_cryst [list $sc_a $sc_b $sc_c 90.0 90.0 90.0]

set ci 0
for {set i 0} {$i < $nx} {incr i} {
    for {set j 0} {$j < $ny} {incr j} {
        for {set k 0} {$k < $nz} {incr k} {
            $base set {x y z} $x0
            $base set resid $r0
            set mb_vec [list [expr $i*$a] [expr $j*$b] [expr $k*$c]]
            $base moveby $mb_vec
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
topology ./ub7.str

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
        }
    }
}
writepsf my_ub7_supercell.psf
writepdb my_ub7_supercell_tmp.pdb

mol new my_ub7_supercell_tmp.pdb
pbc set $sc_cryst
[atomselect top all] writepdb my_ub7_supercell.pdb
exit
