# Full test script for CNDB plugin with trajectory
# Usage: vmd -dispdev text -e test_cndb_full.tcl

puts "Loading CNDB plugin..."
set plugindir [file dirname [info script]]
vmd_plugin_scandirectory $plugindir *.so

puts "\n========================================="
puts "Testing CNDB Plugin with Trajectory"
puts "=========================================\n"

set testfile "/Users/deb0054/github/NDB-Converters/Bcell_short.cndb"

# Load structure and first 10 frames
puts "Loading structure and first 10 frames..."
mol new $testfile type cndb first 0 last 9 waitfor all

set molid [molinfo top]
set natoms [molinfo $molid get numatoms]
set nframes [molinfo $molid get numframes]

puts "\nMolecule loaded successfully!"
puts "  Number of atoms: $natoms"
puts "  Number of frames loaded: $nframes"

# Test coordinates across frames
puts "\nCoordinates of first atom across frames:"
for {set i 0} {$i < $nframes} {incr i} {
    animate goto $i
    set sel [atomselect $molid "index 0"]
    set coords [$sel get {x y z}]
    puts "  Frame $i: $coords"
    $sel delete
}

# Measure distance between two atoms across frames
if {$natoms >= 100} {
    puts "\nDistance between atoms 0 and 99 across frames:"
    for {set i 0} {$i < $nframes} {incr i} {
        animate goto $i
        set sel1 [atomselect $molid "index 0"]
        set sel2 [atomselect $molid "index 99"]
        set pos1 [lindex [$sel1 get {x y z}] 0]
        set pos2 [lindex [$sel2 get {x y z}] 0]
        set dist [veclength [vecsub $pos2 $pos1]]
        puts "  Frame $i: [format %.3f $dist] A"
        $sel1 delete
        $sel2 delete
    }
}

# Test loop constraints
puts "\nChecking loop constraints (bonds)..."
set sel [atomselect $molid "all"]
set bondlist [$sel getbonds]
set total_bonds 0
foreach bonds $bondlist {
    set total_bonds [expr $total_bonds + [llength $bonds]]
}
puts "  Total bonds (chromatin loops): [expr $total_bonds / 2]"
puts "  First 10 atoms and their bonds:"
for {set i 0} {$i < 10} {incr i} {
    set bonds [lindex $bondlist $i]
    if {[llength $bonds] > 0} {
        puts "    Atom $i bonded to: $bonds"
    }
}
$sel delete

puts "\n========================================="
puts "Full trajectory test completed!"
puts "=========================================\n"

# Quit
quit
