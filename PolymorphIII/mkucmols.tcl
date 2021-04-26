package require pbctools

# obabel 154875.cif -opdb -O uc.pdb --fillUC keepconnect
mol new uc.pdb  # contains two molecules
set p [pbc get]
set d2r [expr $M_PI/180.0]
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

set indices1 {0 2 3 5 6 7 8 20 23 24 25 37}
dict set atomnames1 0 "O1"
dict set atomnames1 2 "C6"
dict set atomnames1 3 "C1"
dict set atomnames1 5 "C7"
dict set atomnames1 6 "O3"
dict set atomnames1 7 "O4"
dict set atomnames1 8 "O2"
dict set atomnames1 20 "C3"
dict set atomnames1 23 "C2"
dict set atomnames1 24 "C4"
dict set atomnames1 25 "C5"
dict set atomnames1 37 "C8"
set names [list]
foreach id [dict keys $atomnames1] {
    lappend names [dict get $atomnames1 $id]
}
set mol1 [atomselect top "index $indices1"]
$mol1 set resname "UB7"
$mol1 set resid 1
$mol1 set name $names
set mol1_wrapatoms [atomselect top "index 0 2 3 5 8 25"]
# wrap +a +c
set mv [vecadd $ucellbasisvector1 $ucellbasisvector3]
$mol1_wrapatoms moveby $mv
$mol1 writepdb "unitcell_mol1.pdb"
exit
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
