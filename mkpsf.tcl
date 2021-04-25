package require psfgen
topology /home/cfa/charmm/toppar/top_all36_prot.rtf
topology /home/cfa/charmm/toppar/top_all36_carb.rtf
topology /home/cfa/charmm/toppar/top_all36_lipid.rtf
topology /home/cfa/charmm/toppar/top_all36_na.rtf
topology /home/cfa/charmm/toppar/top_all36_cgenff.rtf
topology ./ub7-acid.str

segment L {
    pdb ub7-acid.pdb
}

coordpdb ub7-acid.pdb L

writepsf my_ub7.psf
writepdb my_ub7.pdb
