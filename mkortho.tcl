package require pbctools
mol new 1x1x1_p.pdb
pbc wrap -orthorhombic
set p [lindex [pbc get] 0]
set c [lindex $p 2]
set beta [lindex $p 4]
set o_c [expr $c*sin($beta*3.141593/180)]
lset p 2 $o_c
lset p 4 90.000
pbc set $p
pbc wrap -shiftcenter {4.5 2.5 2}
set a [atomselect top all]
$a moveby {-4.5 -2.5 -2}
$a set resname UB7

$a writepdb "1x1x1_o.pdb"

exit
