package require pbctools
mol new ortho_unitcell_mol1.pdb
set p [pbc get]
set a [lindex [lindex $p 0] 0]
set b [lindex [lindex $p 0] 1]
set c [lindex [lindex $p 0] 2]

package require psfgen
topology /home/cfa/charmm/toppar/top_all36_prot.rtf
topology /home/cfa/charmm/toppar/top_all36_carb.rtf
topology /home/cfa/charmm/toppar/top_all36_lipid.rtf
topology /home/cfa/charmm/toppar/top_all36_na.rtf
topology /home/cfa/charmm/toppar/top_all36_cgenff.rtf
topology ./ub7.str

segment A {
    pdb ortho_unitcell_mol1.pdb
}
segment B {
    pdb ortho_unitcell_mol2.pdb
}

coordpdb ortho_unitcell_mol1.pdb A
coordpdb ortho_unitcell_mol2.pdb B

guesscoord

writepsf my_ub7_unitcell.psf
writepdb my_ub7_unitcell_tmp.pdb

mol new my_ub7_unitcell_tmp.pdb
pbc set $p
[atomselect top all] writepdb my_ub7_unitcell.pdb
exit
