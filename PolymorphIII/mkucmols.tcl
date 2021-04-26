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

# molecule 1
dict set atomnames1 0 "O1"
dict set atomnames1 1 "H1"
dict set atomnames1 2 "C1"
dict set atomnames1 3 "C2"
dict set atomnames1 4 "H2"
dict set atomnames1 5 "C7"
dict set atomnames1 6 "O3"
dict set atomnames1 7 "O4"
dict set atomnames1 8 "O2"
dict set atomnames1 14 "H4"
dict set atomnames1 20 "C4"
dict set atomnames1 23 "C3"
dict set atomnames1 24 "C5"
dict set atomnames1 25 "C6"
dict set atomnames1 30 "H3"
dict set atomnames1 31 "H5"
dict set atomnames1 32 "H6"
dict set atomnames1 37 "C8"
set indices [list]
set names [list]
foreach id [dict keys $atomnames1] {
    lappend indices $id
    lappend names [dict get $atomnames1 $id]
}
set mol1 [atomselect top "index $indices"]
$mol1 set resname "UB7"
$mol1 set resid 1
$mol1 set name $names
set mol1_wrapatoms [atomselect top "index 0 1 2 3 4 5 8 25 32"]
# wrap +a +c
set mv [vecadd $ucellbasisvector1 $ucellbasisvector3]
$mol1_wrapatoms moveby $mv
$mol1 writepdb "unitcell_mol1.pdb"

# molecule 2
dict set atomnames2 9 "O1"
dict set atomnames2 16 "H1"
dict set atomnames2 38 "C7"
dict set atomnames2 21 "C1"
dict set atomnames2 26 "C2"
dict set atomnames2 33 "H2"
dict set atomnames2 27 "C3"
dict set atomnames2 34 "H3"
dict set atomnames2 22 "C4"
dict set atomnames2 28 "C5"
dict set atomnames2 35 "H5"
dict set atomnames2 29 "C6"
dict set atomnames2 36 "H6"
dict set atomnames2 39 "C8"
dict set atomnames2 10 "O3"
dict set atomnames2 11 "O4"
dict set atomnames2 18 "H4"
dict set atomnames2 12 "O2"
set names [list]
set indices [list]
foreach id [dict keys $atomnames2] {
    lappend indices $id
    lappend names [dict get $atomnames2 $id]
}
set mol2 [atomselect top "index $indices"]
$mol2 set resname "UB7"
$mol2 set resid 2
$mol2 set name $names
set mol2_wrapatoms_c [atomselect top "index 9 12 16 21 26 29 33 36 38"]
# wrap +c
$mol2_wrapatoms_c moveby $ucellbasisvector3
set mol2_wrapatoms_b [atomselect top "index 11 18 12 28 29 35 36"]
$mol2_wrapatoms_b moveby $ucellbasisvector2
$mol2 writepdb "unitcell_mol2.pdb"

exit
