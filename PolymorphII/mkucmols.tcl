package require pbctools

# obabel 1269123.cif -opdb -O uc.pdb --fillUC keepconnect
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
dict set atomnames1 4 "O2"
dict set atomnames1 5 "O1"
dict set atomnames1 6 "C5"
dict set atomnames1 7 "C3"
dict set atomnames1 8 "C4"
dict set atomnames1 9 "C8"
dict set atomnames1 10 "O3"
dict set atomnames1 11 "O4"
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
[atomselect top "index 8 11"] moveby [vecscale $ucellbasisvector1 -1]
[atomselect top "index 0 3"] moveby [vecscale $ucellbasisvector2 -1]
[atomselect top "index 5"] moveby [vecscale [vecadd $ucellbasisvector2 $ucellbasisvector3] -1]
[atomselect top "index 9"] moveby [vecscale [vecadd $ucellbasisvector1 $ucellbasisvector3] -1]
[atomselect top "index 2 7 10"] moveby [vecscale [vecadd $ucellbasisvector1 [vecadd $ucellbasisvector2 $ucellbasisvector3]] -1]

$mol1 writepdb "unitcell_mol1.pdb"
exit
