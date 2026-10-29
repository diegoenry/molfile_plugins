# compare_with_pdb.tcl -- check the pdbx plugin against VMD's PDB reader
#
# For each entry, loads ID.pdb with VMD's pdb plugin and ID.cif with the pdbx
# plugin, then compares atom counts, frame counts and per-atom fields.
#
# Usage (headless):
#   vmd -dispdev text -e tests/compare_with_pdb.tcl -args \
#       /path/to/pdbxplugin.so /path/to/data 1crn 4hhb 1d3z ...
#
# The plugin is loaded with "plugin dlopen"; VMD then uses it in place of the
# built-in pdbx plugin because its version number is higher.

proc field_diffs {a b {tol ""}} {
    set n [llength $a]
    set bad {}
    for {set i 0} {$i < $n} {incr i} {
        set x [lindex $a $i]
        set y [lindex $b $i]
        if {$tol ne ""} {
            if {abs($x - $y) > $tol} { lappend bad $i }
        } elseif {$x ne $y} {
            lappend bad $i
        }
    }
    return $bad
}

proc compare_entry {dir id} {
    set pdb [file join $dir $id.pdb]
    set cif [file join $dir $id.cif]
    set mp [mol new $pdb type pdb waitfor all]
    set mc [mol new $cif type pdbx waitfor all]

    set na [molinfo $mp get numatoms]
    set nb [molinfo $mc get numatoms]
    set fa [molinfo $mp get numframes]
    set fb [molinfo $mc get numframes]
    puts [format "CMP %-5s atoms pdb=%d cif=%d  frames pdb=%d cif=%d" $id $na $nb $fa $fb]
    set ok 1
    if {$na != $nb || $fa != $fb} {
        set ok 0
    } else {
        set sa [atomselect $mp all]
        set sb [atomselect $mc all]
        foreach prop {name resname resid chain insertion altloc element occupancy beta} {
            if {$prop eq "occupancy" || $prop eq "beta"} {
                set d [field_diffs [$sa get $prop] [$sb get $prop] 0.006]
            } else {
                set d [field_diffs [$sa get $prop] [$sb get $prop]]
            }
            if {[llength $d] > 0} {
                set ok 0
                set i [lindex $d 0]
                puts [format "    %-9s %6d differ, e.g. atom %d: pdb='%s' cif='%s'" $prop \
                    [llength $d] $i [lindex [$sa get $prop] $i] [lindex [$sb get $prop] $i]]
            }
        }
        for {set f 0} {$f < $fa} {incr f} {
            $sa frame $f
            $sb frame $f
            set d [field_diffs [$sa get {x y z}] [$sb get {x y z}]]
            set maxd 0.0
            foreach p [$sa get {x y z}] q [$sb get {x y z}] {
                foreach u $p v $q {
                    set e [expr {abs($u - $v)}]
                    if {$e > $maxd} { set maxd $e }
                }
            }
            if {$maxd > 0.0011} {
                set ok 0
                puts [format "    frame %d max coordinate difference %.4f" $f $maxd]
            }
        }
        set bpdb 0; set bcif 0
        foreach l [$sa getbonds] { incr bpdb [llength $l] }
        foreach l [$sb getbonds] { incr bcif [llength $l] }
        puts [format "    bonds pdb=%d cif=%d   cell pdb={%s} cif={%s}" \
            [expr {$bpdb / 2}] [expr {$bcif / 2}] \
            [molinfo $mp get {a b c}] [molinfo $mc get {a b c}]]
        $sa delete
        $sb delete
    }
    puts [format "CMP %-5s %s" $id [expr {$ok ? "MATCH" : "DIFF"}]]
    mol delete $mp
    mol delete $mc
}

set plugin_so [lindex $argv 0]
set data_dir  [lindex $argv 1]
plugin dlopen $plugin_so
foreach id [lrange $argv 2 end] {
    if {[catch {compare_entry $data_dir $id} err]} {
        puts "CMP $id ERROR $err"
    }
}
quit
