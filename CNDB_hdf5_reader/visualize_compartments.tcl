# Visualization script for chromatin compartments
# Usage: vmd -e visualize_compartments.tcl -args <cndb_file> [frame]

# Load plugin
set scriptdir [file dirname [info script]]
if {[file exists "$scriptdir/cndbplugin.so"]} {
    vmd_plugin_scandirectory $scriptdir *.so
}

# Get arguments
if {$argc >= 1} {
    set cndb_file [lindex $argv 0]
} else {
    puts "Usage: vmd -e visualize_compartments.tcl -args <cndb_file> \[frame\]"
    puts "Example: vmd -e visualize_compartments.tcl -args data.cndb 0"
    exit
}

set frame_num 0
if {$argc >= 2} {
    set frame_num [lindex $argv 1]
}

puts "=========================================="
puts "Chromatin Compartment Visualization"
puts "=========================================="
puts "File: $cndb_file"
puts "Frame: $frame_num"
puts ""

# Load single frame
mol new $cndb_file type cndb first $frame_num last $frame_num
set molid [molinfo top]

# Get statistics
set natoms [molinfo $molid get numatoms]
puts "Total atoms: $natoms"

# Count each compartment type
foreach type {A1 A2 B1 B3 NA} {
    set sel [atomselect $molid "type $type"]
    set count [$sel num]
    $sel delete
    puts "  $type: $count beads"
}
puts ""

# Delete default representation
mol delrep 0 $molid

# Setup color scheme for compartment types (if in GUI mode)
catch {color Type A1 red}
catch {color Type A2 orange}
catch {color Type B1 blue}
catch {color Type B3 cyan}
catch {color Type NA gray}

# Representation 1: Beads colored by compartment type
mol representation VDW 5.0 12.0
mol color Type
mol selection "all"
mol material Opaque
mol addrep $molid
puts "Created representation 1: All beads colored by compartment type"

# Representation 2: Chromatin loops
mol representation Bonds 2.0 12.0
mol color ColorID 10
mol selection "all"
mol material Opaque
mol addrep $molid
puts "Created representation 2: Chromatin loops (tan)"

# Representation 3: Active chromatin only (A compartments)
mol representation VDW 6.0 12.0
mol color Type
mol selection "type A1 A2"
mol material Opaque
mol addrep $molid
mol showrep $molid 2 off
puts "Created representation 3: Active chromatin (A1, A2) - OFF by default"

# Representation 4: Inactive chromatin only (B compartments)
mol representation VDW 6.0 12.0
mol color Type
mol selection "type B1 B3"
mol material Opaque
mol addrep $molid
mol showrep $molid 3 off
puts "Created representation 4: Inactive chromatin (B1, B3) - OFF by default"

# Display settings
display projection orthographic
display depthcue off
color Display Background white
axes location off
display resetview

puts ""
puts "=========================================="
puts "Color Scheme:"
puts "  A1 (active):   RED"
puts "  A2 (active):   ORANGE"
puts "  B1 (inactive): BLUE"
puts "  B3 (inactive): CYAN"
puts "  NA:            GRAY"
puts "=========================================="
puts ""
puts "Visualization tips:"
puts "  - Representations 3 & 4 show only A or B compartments"
puts "  - Toggle them in Graphics > Representations"
puts "  - Use 'mol on 2' or 'mol off 2' to toggle rep 3"
puts "=========================================="
