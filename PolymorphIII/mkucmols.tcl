package require pbctools

mol new 1x1x1.pdb  # contains two molecules
set p [pbc get]

set heavies [atomselect top "noh"]

set mol1 [atomselect top "index 0 1 2 3 16 17 20 21 22 23 36 37"]
$mol1 set resname "UB7"
$mol1 set resid 1
$mol1 set name {"O1" "O3" "O4" "O2" "C6" "C3" "C1" "C2" "C4" "C5" "C7" "C8"}

set mol2 [atomselect top "index 4 5 6 7 18 19 24 25 26 27 38 39"]
$mol2 set resname "UB7"
$mol2 set resid 2
$mol2 set name {"O1" "O3" "O4" "O2" "C6" "C3" "C1" "C2" "C4" "C5" "C7" "C8"}

set lx [lindex [lindex $p 0] 0]
set moveus [atomselect top "resid 2 and name C1 C5 C6 C7 O1 O2"]
$moveus moveby [list $lx 0 0]

$mol1 writepdb "unitcell_mol1.pdb"
$mol2 writepdb "unitcell_mol2.pdb"

exit
