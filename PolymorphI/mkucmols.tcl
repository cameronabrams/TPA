package require pbctools

# obabel 1269122.cif -opdb -O uc.pdb --fillUC keepconnect
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
dict set atomnames1 0 "C1"
dict set atomnames1 1 "C6"
dict set atomnames1 2 "C2"
dict set atomnames1 3 "C7"
dict set atomnames1 4 "H1"
dict set atomnames1 5 "H6"
dict set atomnames1 6 "H2"
dict set atomnames1 7 "O1"
dict set atomnames1 8 "O2"
dict set atomnames1 9 "C5"
dict set atomnames1 10 "C3"
dict set atomnames1 11 "C4"
dict set atomnames1 12 "H5"
dict set atomnames1 13 "H3"
dict set atomnames1 14 "C8"
dict set atomnames1 15 "O4"
dict set atomnames1 16 "O3"
dict set atomnames1 17 "H4"
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
# wrap +b
[atomselect top "index 0 3 7 4"] moveby $ucellbasisvector2
[atomselect top "index 11 14 15 17"] moveby [vecscale [vecadd $ucellbasisvector1 $ucellbasisvector3] -1]
[atomselect top "index 1 5 9 12"] moveby [vecadd $ucellbasisvector2 [vecscale [vecadd $ucellbasisvector1 $ucellbasisvector3] -1]]
[atomselect top "index 8"] moveby [vecadd $ucellbasisvector2 [vecscale $ucellbasisvector3 -1]]
[atomselect top "index 16"] moveby [vecscale $ucellbasisvector1 -1]

$mol1 writepdb "unitcell_mol1.pdb"
exit
