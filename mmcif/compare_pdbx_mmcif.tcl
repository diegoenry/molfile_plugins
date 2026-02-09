# compare_pdbx_mmcif.tcl
#
# Compare molecules loaded with VMD's built-in pdbx plugin vs. the custom
# mmcif molfile plugin.  Reports differences in atom properties, bonds,
# secondary structure, coordinates, and unit-cell parameters.
#
# Usage inside VMD:
#   cd /path/to/mmcif
#   source compare_pdbx_mmcif.tcl
#   compare_plugins 8EUU-assembly1.cif
#
# The mmcifplugin.so must be in the current directory.

# ── Helpers ──────────────────────────────────────────────────────────────

proc section {title} {
    puts "\n[string repeat = 60]"
    puts "  $title"
    puts [string repeat = 60]
}

proc subsection {title} {
    puts "\n--- $title ---"
}

# Return a dict summarising one loaded molecule.
# Keys: molid, natoms, plus per-atom lists for every property we compare.
proc collect_mol_data {molid} {
    set d [dict create molid $molid]
    dict set d natoms [molinfo $molid get numatoms]

    # Atom-level properties available in VMD
    set sel [atomselect $molid "all"]
    foreach prop {name type resname resid chain segname altloc insertion
                  occupancy beta mass charge atomicnumber element
                  x y z index} {
        dict set d $prop [$sel get $prop]
    }
    $sel delete

    # Bond topology
    set sel [atomselect $molid "all"]
    dict set d bondlist  [$sel getbonds]
    dict set d bondorder [$sel getbondorders]
    $sel delete

    # Secondary structure
    set sel [atomselect $molid "all"]
    dict set d structure [$sel get structure]
    $sel delete

    # Unit cell
    foreach param {a b c alpha beta gamma} {
        dict set d cell_$param [molinfo $molid get $param]
    }

    # Number of frames (models)
    dict set d nframes [molinfo $molid get numframes]

    return $d
}

# Compare two flat lists element-by-element.  Return indices that differ.
proc compare_lists {listA listB {tol 0.0}} {
    set diffs {}
    set n [llength $listA]
    for {set i 0} {$i < $n} {incr i} {
        set a [lindex $listA $i]
        set b [lindex $listB $i]
        if {$tol > 0.0} {
            if {abs($a - $b) > $tol} { lappend diffs $i }
        } else {
            if {$a ne $b} { lappend diffs $i }
        }
    }
    return $diffs
}

# Pretty-print a handful of differing atoms for a given property.
proc report_diffs {prop diffs dataA dataB {max_show 20}} {
    set nd [llength $diffs]
    if {$nd == 0} {
        puts "  $prop: IDENTICAL"
        return
    }
    puts "  $prop: $nd differences"
    set namesA [dict get $dataA name]
    set residsA [dict get $dataA resid]
    set chainsA [dict get $dataA chain]
    set valsA [dict get $dataA $prop]
    set valsB [dict get $dataB $prop]
    set shown 0
    foreach idx $diffs {
        if {$shown >= $max_show} {
            puts "    ... and [expr {$nd - $max_show}] more"
            break
        }
        set label "[lindex $chainsA $idx]:[lindex $residsA $idx]:[lindex $namesA $idx]"
        puts [format "    atom %5d %-20s  pdbx=%-12s  mmcif=%-12s" \
                  $idx $label [lindex $valsA $idx] [lindex $valsB $idx]]
        incr shown
    }
}

# ── Main comparison procedure ────────────────────────────────────────────

proc compare_plugins {ciffile} {

    # ── 0. Register mmcif plugin ─────────────────────────────────────────
    section "Loading plugins"

    set plugindir [pwd]
    if {[file exists [file join $plugindir mmcifplugin.so]]} {
        vmd_plugin_scandirectory [pwd] mmcifplugin.so 
        puts "  Scanned $plugindir for plugins"
    } else {
        puts "  WARNING: mmcifplugin.so not found in $plugindir"
    }

    # ── 1. Load with pdbx (VMD built-in) ─────────────────────────────────
    section "Loading with pdbx (VMD built-in)"
    mol new $ciffile type pdbx waitfor all
    set mol_pdbx [molinfo top]
    puts "  molid $mol_pdbx — [molinfo $mol_pdbx get numatoms] atoms, [molinfo $mol_pdbx get numframes] frames"

    # ── 2. Load with mmcif (custom plugin) ───────────────────────────────
    section "Loading with mmcif (custom plugin)"
    mol new $ciffile type mmcif waitfor all
    set mol_mmcif [molinfo top]
    puts "  molid $mol_mmcif — [molinfo $mol_mmcif get numatoms] atoms, [molinfo $mol_mmcif get numframes] frames"

    # ── 3. Collect data ──────────────────────────────────────────────────
    section "Collecting atom data"
    set A [collect_mol_data $mol_pdbx]
    set B [collect_mol_data $mol_mmcif]
    puts "  pdbx:  [dict get $A natoms] atoms, [dict get $A nframes] frames"
    puts "  mmcif: [dict get $B natoms] atoms, [dict get $B nframes] frames"

    # ── 4. Atom count ────────────────────────────────────────────────────
    section "Atom Count"
    set na [dict get $A natoms]
    set nb [dict get $B natoms]
    if {$na != $nb} {
        puts "  DIFFERENT: pdbx=$na  mmcif=$nb"
        puts "  Cannot do per-atom comparison with different atom counts."
        puts "  Cleaning up..."
        mol delete $mol_pdbx
        mol delete $mol_mmcif
        return
    }
    puts "  Both have $na atoms — OK"

    # ── 5. Atom identity properties (string comparison) ──────────────────
    section "Atom Identity Properties"
    foreach prop {name type resname element chain segname altloc insertion} {
        set diffs [compare_lists [dict get $A $prop] [dict get $B $prop]]
        report_diffs $prop $diffs $A $B
    }

    # ── 6. Residue IDs ──────────────────────────────────────────────────
    section "Residue IDs"
    set diffs [compare_lists [dict get $A resid] [dict get $B resid]]
    report_diffs resid $diffs $A $B

    # ── 7. Numeric atom properties (with tolerance) ─────────────────────
    section "Numeric Atom Properties"
    foreach {prop tol} {occupancy 0.001 beta 0.01 mass 0.01 charge 0.01} {
        set diffs [compare_lists [dict get $A $prop] [dict get $B $prop] $tol]
        report_diffs $prop $diffs $A $B
    }

    subsection "Atomic numbers"
    set diffs [compare_lists [dict get $A atomicnumber] [dict get $B atomicnumber]]
    report_diffs atomicnumber $diffs $A $B

    # ── 8. Coordinates ──────────────────────────────────────────────────
    section "Coordinates (frame 0)"
    set coord_tol 0.001
    foreach axis {x y z} {
        set diffs [compare_lists [dict get $A $axis] [dict get $B $axis] $coord_tol]
        set nd [llength $diffs]
        if {$nd == 0} {
            puts "  $axis: IDENTICAL (tol=$coord_tol)"
        } else {
            puts "  $axis: $nd atoms differ (tol=$coord_tol)"
            set shown 0
            foreach idx [lrange $diffs 0 4] {
                set va [lindex [dict get $A $axis] $idx]
                set vb [lindex [dict get $B $axis] $idx]
                puts [format "    atom %5d  pdbx=%12.4f  mmcif=%12.4f  delta=%g" \
                          $idx $va $vb [expr {$va - $vb}]]
            }
            if {$nd > 5} { puts "    ... and [expr {$nd - 5}] more" }
        }
    }

    # Compute overall RMSD between the two coordinate sets
    set sum2 0.0
    set xA [dict get $A x]; set yA [dict get $A y]; set zA [dict get $A z]
    set xB [dict get $B x]; set yB [dict get $B y]; set zB [dict get $B z]
    for {set i 0} {$i < $na} {incr i} {
        set dx [expr {[lindex $xA $i] - [lindex $xB $i]}]
        set dy [expr {[lindex $yA $i] - [lindex $yB $i]}]
        set dz [expr {[lindex $zA $i] - [lindex $zB $i]}]
        set sum2 [expr {$sum2 + $dx*$dx + $dy*$dy + $dz*$dz}]
    }
    set rmsd [expr {sqrt($sum2 / $na)}]
    puts [format "\n  Coordinate RMSD between plugins: %.6f A" $rmsd]

    # ── 9. Unit cell ────────────────────────────────────────────────────
    section "Unit Cell Parameters"
    set cell_tol 0.01
    foreach param {a b c alpha beta gamma} {
        set va [dict get $A cell_$param]
        set vb [dict get $B cell_$param]
        set match [expr {abs($va - $vb) < $cell_tol}]
        if {$match} {
            puts [format "  %-6s: %10.4f  (both agree)" $param $va]
        } else {
            puts [format "  %-6s: pdbx=%10.4f  mmcif=%10.4f  DIFFERENT" $param $va $vb]
        }
    }

    # ── 10. Secondary structure ─────────────────────────────────────────
    section "Secondary Structure"
    set ssA [dict get $A structure]
    set ssB [dict get $B structure]

    # Count by type for each plugin
    foreach lbl {pdbx mmcif} ss [list $ssA $ssB] {
        set counts [dict create]
        foreach s $ss {
            dict incr counts $s
        }
        puts "  $lbl SS distribution:"
        foreach key [lsort [dict keys $counts]] {
            set label $key
            switch $key {
                C { set label "C (coil)" }
                H { set label "H (helix)" }
                E { set label "E (strand)" }
                T { set label "T (turn)" }
                G { set label "G (3-10 helix)" }
                I { set label "I (pi helix)" }
                B { set label "B (bridge)" }
            }
            puts [format "    %-16s %6d atoms" $label [dict get $counts $key]]
        }
    }

    set diffs [compare_lists $ssA $ssB]
    set nd [llength $diffs]
    if {$nd == 0} {
        puts "\n  Secondary structure: IDENTICAL"
    } else {
        puts "\n  Secondary structure: $nd differences out of $na atoms"
        # Summarise: how many switched from what to what
        set transitions [dict create]
        foreach idx $diffs {
            set key "[lindex $ssA $idx]->[lindex $ssB $idx]"
            dict incr transitions $key
        }
        puts "  Transition summary (pdbx -> mmcif):"
        foreach key [lsort [dict keys $transitions]] {
            puts [format "    %-10s %5d atoms" $key [dict get $transitions $key]]
        }
        # Show a few examples
        puts "  Examples:"
        set shown 0
        set namesA [dict get $A name]
        set residsA [dict get $A resid]
        set chainsA [dict get $A chain]
        set resnamesA [dict get $A resname]
        foreach idx $diffs {
            if {$shown >= 10} {
                puts "    ... and [expr {$nd - 10}] more"
                break
            }
            set label "[lindex $chainsA $idx]:[lindex $resnamesA $idx][lindex $residsA $idx]:[lindex $namesA $idx]"
            puts [format "    atom %5d %-20s  pdbx=%-4s  mmcif=%-4s" \
                      $idx $label [lindex $ssA $idx] [lindex $ssB $idx]]
            incr shown
        }
    }

    # ── 11. Bonds ───────────────────────────────────────────────────────
    section "Bond Topology"

    set bondsA [dict get $A bondlist]
    set bondsB [dict get $B bondlist]
    set ordersA [dict get $A bondorder]
    set ordersB [dict get $B bondorder]

    # Count total bonds (each bond counted once per atom that lists it)
    set nbA 0; foreach bl $bondsA { incr nbA [llength $bl] }
    set nbB 0; foreach bl $bondsB { incr nbB [llength $bl] }
    puts "  pdbx  total bond entries (sum of per-atom lists): $nbA  ([expr {$nbA/2}] unique bonds)"
    puts "  mmcif total bond entries (sum of per-atom lists): $nbB  ([expr {$nbB/2}] unique bonds)"

    # Build canonical bond sets for comparison: set of "i-j" with i<j
    proc build_bond_set {bondlists} {
        set bset [dict create]
        set aidx 0
        foreach bl $bondlists {
            foreach partner $bl {
                set i $aidx
                set j $partner
                if {$i > $j} { set i $partner; set j $aidx }
                dict set bset "$i-$j" 1
            }
            incr aidx
        }
        return $bset
    }

    proc build_bond_order_map {bondlists orderlists} {
        set bmap [dict create]
        set aidx 0
        foreach bl $bondlists ol $orderlists {
            set k 0
            foreach partner $bl {
                set i $aidx
                set j $partner
                if {$i > $j} { set i $partner; set j $aidx }
                set order [lindex $ol $k]
                dict set bmap "$i-$j" $order
                incr k
            }
            incr aidx
        }
        return $bmap
    }

    set setA [build_bond_set $bondsA]
    set setB [build_bond_set $bondsB]
    set nA [dict size $setA]
    set nB [dict size $setB]
    puts "\n  Unique bonds: pdbx=$nA  mmcif=$nB"

    # Bonds in pdbx but not mmcif
    set only_pdbx {}
    dict for {bond _} $setA {
        if {![dict exists $setB $bond]} { lappend only_pdbx $bond }
    }
    # Bonds in mmcif but not pdbx
    set only_mmcif {}
    dict for {bond _} $setB {
        if {![dict exists $setA $bond]} { lappend only_mmcif $bond }
    }
    # Common bonds
    set common {}
    dict for {bond _} $setA {
        if {[dict exists $setB $bond]} { lappend common $bond }
    }

    puts "  Common bonds:     [llength $common]"
    puts "  Only in pdbx:     [llength $only_pdbx]"
    puts "  Only in mmcif:    [llength $only_mmcif]"

    # Print a few examples of bonds unique to each
    set namesA [dict get $A name]
    set residsA [dict get $A resid]
    set chainsA [dict get $A chain]
    set resnamesA [dict get $A resname]

    proc atom_label {names resnames resids chains idx} {
        return "[lindex $chains $idx]:[lindex $resnames $idx][lindex $resids $idx]:[lindex $names $idx]"
    }

    if {[llength $only_pdbx] > 0} {
        subsection "Bonds only in pdbx (first 20)"
        set shown 0
        foreach bond [lrange $only_pdbx 0 19] {
            set parts [split $bond -]
            set i [lindex $parts 0]; set j [lindex $parts 1]
            set li [atom_label $namesA $resnamesA $residsA $chainsA $i]
            set lj [atom_label $namesA $resnamesA $residsA $chainsA $j]
            puts [format "    %s -- %s  (atoms %d-%d)" $li $lj $i $j]
            incr shown
        }
        if {[llength $only_pdbx] > 20} {
            puts "    ... and [expr {[llength $only_pdbx] - 20}] more"
        }
    }

    if {[llength $only_mmcif] > 0} {
        subsection "Bonds only in mmcif (first 20)"
        set shown 0
        foreach bond [lrange $only_mmcif 0 19] {
            set parts [split $bond -]
            set i [lindex $parts 0]; set j [lindex $parts 1]
            set li [atom_label $namesA $resnamesA $residsA $chainsA $i]
            set lj [atom_label $namesA $resnamesA $residsA $chainsA $j]
            puts [format "    %s -- %s  (atoms %d-%d)" $li $lj $i $j]
        }
        if {[llength $only_mmcif] > 20} {
            puts "    ... and [expr {[llength $only_mmcif] - 20}] more"
        }
    }

    # Compare bond orders for common bonds
    if {[llength $common] > 0} {
        subsection "Bond Order Comparison (common bonds)"
        set mapA [build_bond_order_map $bondsA $ordersA]
        set mapB [build_bond_order_map $bondsB $ordersB]
        set order_diffs {}
        foreach bond $common {
            set oA [dict get $mapA $bond]
            set oB [dict get $mapB $bond]
            # Treat -1 (unknown) specially
            if {$oA ne $oB} {
                lappend order_diffs [list $bond $oA $oB]
            }
        }
        if {[llength $order_diffs] == 0} {
            puts "  Bond orders: IDENTICAL for all [llength $common] common bonds"
        } else {
            puts "  Bond orders differ for [llength $order_diffs] common bonds"
            set shown 0
            foreach diff [lrange $order_diffs 0 19] {
                lassign $diff bond oA oB
                set parts [split $bond -]
                set i [lindex $parts 0]; set j [lindex $parts 1]
                set li [atom_label $namesA $resnamesA $residsA $chainsA $i]
                set lj [atom_label $namesA $resnamesA $residsA $chainsA $j]
                puts [format "    %s -- %s  pdbx_order=%-6s mmcif_order=%-6s" $li $lj $oA $oB]
            }
            if {[llength $order_diffs] > 20} {
                puts "    ... and [expr {[llength $order_diffs] - 20}] more"
            }
        }
    }

    # ── 12. Per-residue bond count comparison ───────────────────────────
    subsection "Per-atom bond count differences"
    set bond_count_diffs {}
    for {set i 0} {$i < $na} {incr i} {
        set cntA [llength [lindex $bondsA $i]]
        set cntB [llength [lindex $bondsB $i]]
        if {$cntA != $cntB} {
            lappend bond_count_diffs [list $i $cntA $cntB]
        }
    }
    if {[llength $bond_count_diffs] == 0} {
        puts "  Per-atom bond counts: IDENTICAL"
    } else {
        puts "  [llength $bond_count_diffs] atoms have different bond counts"
        foreach diff [lrange $bond_count_diffs 0 19] {
            lassign $diff idx cntA cntB
            set label [atom_label $namesA $resnamesA $residsA $chainsA $idx]
            puts [format "    atom %5d %-20s  pdbx=%d bonds  mmcif=%d bonds" \
                      $idx $label $cntA $cntB]
        }
        if {[llength $bond_count_diffs] > 20} {
            puts "    ... and [expr {[llength $bond_count_diffs] - 20}] more"
        }
    }

    # ── 13. Summary ─────────────────────────────────────────────────────
    section "SUMMARY"
    puts "  File: $ciffile"
    puts "  Atoms: $na (both plugins agree)"
    puts [format "  Coordinate RMSD: %.6f A" $rmsd]
    puts "  Frames: pdbx=[dict get $A nframes]  mmcif=[dict get $B nframes]"

    set total_prop_diffs 0
    foreach prop {name type resname resid element chain segname altloc insertion} {
        set nd [llength [compare_lists [dict get $A $prop] [dict get $B $prop]]]
        if {$nd > 0} {
            puts "  $prop: $nd differences"
            incr total_prop_diffs $nd
        }
    }
    if {$total_prop_diffs == 0} {
        puts "  Atom properties: ALL IDENTICAL"
    }

    set ss_diffs [llength [compare_lists $ssA $ssB]]
    puts "  Secondary structure differences: $ss_diffs"
    puts "  Bonds only in pdbx:  [llength $only_pdbx]"
    puts "  Bonds only in mmcif: [llength $only_mmcif]"
    puts ""

    # Clean up
    mol delete $mol_pdbx
    mol delete $mol_mmcif
}

# ── Auto-run if a default test file is present ──────────────────────────
puts "compare_pdbx_mmcif.tcl loaded."
puts "Usage: compare_plugins <file.cif>"
puts ""

if {[file exists 8EUU-assembly1.cif]} {
    puts "Found 8EUU-assembly1.cif — running comparison automatically...\n"
    compare_plugins 8EUU-assembly1.cif
}
