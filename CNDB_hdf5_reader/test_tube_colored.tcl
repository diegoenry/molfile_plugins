# Test script for chromatin tube with color-by-type

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

puts "\n=========================================="
puts "Test: Tube colored by compartment type"
puts "=========================================="

# Draw chromatin tube colored by type (using defaults for radius, segments, resolution)
draw_chromatin_tube $molid 10.0 12 6 blue Opaque 1

puts "\n=========================================="
puts "Test complete!"
puts "=========================================="
puts "\nYou should see a smooth tube with colors:"
puts "  RED    - A1 compartments"
puts "  ORANGE - A2 compartments"
puts "  BLUE   - B1 compartments"
puts "  CYAN   - B3 compartments"
puts "  GRAY   - NA regions"
puts ""
puts "To clear: graphics top delete all"
puts ""
