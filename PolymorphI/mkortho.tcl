package require pbctools
mol new 1x1x1_p.pdb
[atomselect top "index 0 to 17"] writepdb 1x1x1_p1.pdb
mol delete top
mol new 1x1x1_p1.pdb

pbc wrap -orthorhombic
set p [lindex [pbc get] 0]
set a [lindex $p 0]
set b [lindex $p 1]
set c [lindex $p 2]
set alpha [lindex $p 3]
set beta [lindex $p 4]
set gamma [lindex $p 5]
set o_a [expr $a*sin($gamma*3.141593/180)]
set o_c [expr $c*sin($beta*3.141593/180)*sin($alpha*3.141593/180)]
lset p 0 $o_a
lset p 2 $o_c
lset p 3 90.000
lset p 4 90.000
lset p 5 90.000
pbc set $p
pbc wrap -shiftcenter {4.5 2.5 2}
set a [atomselect top all]
$a moveby {-4.5 -2.5 -2}
$a set resname UB7

$a writepdb "1x1x1_o.pdb"

exit
