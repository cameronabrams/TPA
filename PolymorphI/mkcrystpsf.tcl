package require pbctools
mol new unitcell_mol1.pdb
set p [pbc get]

package require psfgen
topology /home/cfa/charmm/toppar/top_all36_prot.rtf
topology /home/cfa/charmm/toppar/top_all36_carb.rtf
topology /home/cfa/charmm/toppar/top_all36_lipid.rtf
topology /home/cfa/charmm/toppar/top_all36_na.rtf
topology /home/cfa/charmm/toppar/top_all36_cgenff.rtf
topology ../charmm/ub7.str

segment A {
    pdb unitcell_mol1.pdb
}

coordpdb unitcell_mol1.pdb A

guesscoord

writepsf my_ub7_unitcell.psf
writepdb my_ub7_unitcell_tmp.pdb

mol new my_ub7_unitcell_tmp.pdb
pbc set $p
[atomselect top all] writepdb my_ub7_unitcell.pdb
exit
