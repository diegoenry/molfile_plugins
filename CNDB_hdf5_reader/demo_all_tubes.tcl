# Demonstration of all tube drawing modes
# Shows the three different ways to visualize chromatin tubes

# Load plugin and tube functions
set scriptdir [file dirname [info script]]
if {[file exists "$scriptdir/cndbplugin.so"]} {
    vmd_plugin_scandirectory $scriptdir *.so
}
source "$scriptdir/chromatin_tube.tcl"

# Load CNDB file
set testfile "/Users/deb0054/github/NDB-Converters/Bcell_short.cndb"

puts "\n=========================================="
puts "Chromatin Tube Visualization Demo"
puts "=========================================="
puts "This demo shows three visualization modes:"
puts "  1. Single color tube"
puts "  2. Rainbow tube (colored by type)"
puts "  3. Separate tubes per type"
puts "=========================================="

# Mode 1: Single color tube
puts "\nMode 1: Single Color Tube"
puts "------------------------------------------"
mol new $testfile type cndb first 0 last 0
set mol1 [molinfo top]
mol rename $mol1 "Single Color"
mol delrep 0 $mol1

draw_chromatin_tube $mol1

puts "\nMolecule 0: Single blue tube created"

# Mode 2: Rainbow tube (colored by type)
puts "\n\nMode 2: Rainbow Tube (Color by Type)"
puts "------------------------------------------"
mol new $testfile type cndb first 0 last 0
set mol2 [molinfo top]
mol rename $mol2 "Rainbow Tube"
mol delrep 0 $mol2

draw_chromatin_tube $mol2 10.0 12 6 blue Opaque 1

puts "\nMolecule 1: Rainbow tube created (colored by compartment)"

# Mode 3: Separate tubes per type
puts "\n\nMode 3: Separate Tubes per Type"
puts "------------------------------------------"
mol new $testfile type cndb first 0 last 0
set mol3 [molinfo top]
mol rename $mol3 "Separate Tubes"
mol delrep 0 $mol3

draw_chromatin_tube_by_type $mol3

puts "\nMolecule 2: Separate tubes for each compartment type"

puts "\n=========================================="
puts "Demo Complete!"
puts "=========================================="
puts "\nThree molecules loaded:"
puts "  Mol 0: Single color tube (blue)"
puts "  Mol 1: Rainbow tube (A1=red, A2=orange, B1=blue, B3=cyan, NA=gray)"
puts "  Mol 2: Separate tubes per compartment type"
puts ""
puts "Toggle visibility in VMD:"
puts "  mol on/off 0   - Show/hide single color"
puts "  mol on/off 1   - Show/hide rainbow tube"
puts "  mol on/off 2   - Show/hide separate tubes"
puts ""
puts "Clear all graphics:"
puts "  graphics 0 delete all"
puts "  graphics 1 delete all"
puts "  graphics 2 delete all"
puts "=========================================="
