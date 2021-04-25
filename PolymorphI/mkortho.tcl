package require pbctools
mol new 1x1x1_p.pdb
[atomselect top "index 0 to 17"] writepdb 1x1x1_p1.pdb
mol delete top
mol new 1x1x1_p1.pdb

set p [lindex [pbc get] 0]
set a [lindex $p 0]
set b [lindex $p 1]
set c [lindex $p 2]
set d2r [expr 3.141593/180.0]
set alpha [expr [lindex $p 3]*$d2r]
set beta [expr [lindex $p 4]*$d2r]
set gamma [expr [lindex $p 5]*$d2r]
set o_a $a
set xy [expr $b*cos($gamma)]
set xz [expr $c*cos($beta)]
set ly [expr sqrt($b*$b-$xy*$xy)]
set yz [expr ($b*$c*cos($alpha)-$xy*$xz)/$ly]
set lz [expr sqrt($c*$c-$xz*$xz-$yz*$yz)]
set o_b $ly
set o_c $lz
lset p 0 $o_a
lset p 1 $o_b
lset p 2 $o_c
lset p 3 90.000
lset p 4 90.000
lset p 5 90.000
pbc set $p
#pbc wrap
#pbc wrap -shiftcenter {4.5 2.5 2}
set a [atomselect top all]
#$a moveby {-4.5 -2.5 -2}
#$a set resname UB7

$a writepdb "1x1x1_o.pdb"

exit
