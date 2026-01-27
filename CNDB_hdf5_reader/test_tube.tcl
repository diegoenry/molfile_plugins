# Test script for chromatin tube visualization

# Load plugin
set scriptdir [file dirname [info script]]
if {[file exists "$scriptdir/cndbplugin.so"]} {
    vmd_plugin_scandirectory $scriptdir *.so
}

# Load chromatin tube functions
source "$scriptdir/chromatin_tube.tcl"

# Load CNDB file (single frame)
set testfile "/Users/deb0054/github/NDB-Converters/Bcell_short.cndb"
mol new $testfile type cndb first 0 last 0

set molid [molinfo top]

# Delete default representation
mol delrep 0 $molid

# Draw chromatin tube
puts "\n=========================================="
puts "Test 1: Single tube through all beads"
puts "=========================================="
draw_chromatin_tube $molid 10.0 16 6 blue Opaque

# Wait and clear
puts "\nPress Enter to continue to test 2..."
gets stdin

graphics $molid delete all

# Draw tubes by compartment type
puts "\n=========================================="
puts "Test 2: Tubes colored by compartment type"
puts "=========================================="
draw_chromatin_tube_by_type $molid 8.0 12 4 Opaque

puts "\n=========================================="
puts "Tests complete!"
puts "=========================================="
puts "\nTo manually clear graphics:"
puts "  graphics top delete all"
puts ""

# Don't quit automatically - let user examine
