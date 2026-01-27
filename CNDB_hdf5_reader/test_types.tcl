# Test script to check bead types
set scriptdir [file dirname [info script]]
if {[file exists "$scriptdir/cndbplugin.so"]} {
    vmd_plugin_scandirectory $scriptdir *.so
}

set testfile "/Users/deb0054/github/NDB-Converters/Bcell_short.cndb"
mol new $testfile type cndb first 0 last 0

set molid [molinfo top]

# Get all unique types
set sel [atomselect $molid "all"]
set types [$sel get type]
set names [$sel get name]
$sel delete

# Count occurrences of each type
array set type_counts {}
foreach type $types {
    if {[info exists type_counts($type)]} {
        incr type_counts($type)
    } else {
        set type_counts($type) 1
    }
}

puts "\n=========================================="
puts "Bead Type Distribution"
puts "=========================================="
foreach type [lsort [array names type_counts]] {
    puts "  $type: $type_counts($type) beads"
}
puts "=========================================="

# Show first 20 atoms
puts "\nFirst 20 atoms:"
puts "  Index  Name  Type  Resname  Genomic_Start"
for {set i 0} {$i < 20} {incr i} {
    set sel [atomselect $molid "index $i"]
    set name [lindex [$sel get name] 0]
    set type [lindex [$sel get type] 0]
    set resname [lindex [$sel get resname] 0]
    set beta [lindex [$sel get beta] 0]
    puts [format "  %5d  %-4s  %-4s  %-7s  %.0f" $i $name $type $resname $beta]
    $sel delete
}

quit
