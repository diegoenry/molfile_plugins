# Script to create publication-quality figures from CNDB files
# Usage: vmd -dispdev text -e make_figure.tcl -args input.cndb output_basename [frame]

# Load plugin
set scriptdir [file dirname [info script]]
if {[file exists "$scriptdir/cndbplugin.so"]} {
    vmd_plugin_scandirectory $scriptdir *.so
}

# Get arguments
if {$argc < 2} {
    puts "Usage: vmd -dispdev text -e make_figure.tcl -args <input.cndb> <output_basename> \[frame\]"
    puts "Example: vmd -dispdev text -e make_figure.tcl -args data.cndb chromatin_structure 0"
    exit
}

set input_file [lindex $argv 0]
set output_base [lindex $argv 1]
set frame_num 0
if {$argc >= 3} {
    set frame_num [lindex $argv 2]
}

puts "Creating figure from CNDB file..."
puts "  Input: $input_file"
puts "  Output: ${output_base}_*.tga"
puts "  Frame: $frame_num"

# Load molecule (only one frame)
mol new $input_file type cndb first $frame_num last $frame_num

set molid [molinfo top]

# Delete default representation
mol delrep 0 $molid

# Setup 1: Beads colored by genomic position with loops
mol representation VDW 5.0 12.0
mol color Beta
mol selection "all"
mol material Opaque
mol addrep $molid

mol representation Bonds 2.0 12.0
mol color ColorID 3
mol selection "all"
mol material Opaque
mol addrep $molid

# Display settings
display projection orthographic
display depthcue off
color Display Background white
axes location off
display resetview

# Render view 1
puts "Rendering view 1 (genomic position coloring + loops)..."
render TachyonInternal ${output_base}_view1.tga

# Change to different coloring scheme
mol modcolor 0 $molid Occupancy
puts "Rendering view 2 (occupancy coloring)..."
render TachyonInternal ${output_base}_view2.tga

# Setup 2: Just beads, larger
mol delete $molid
mol new $input_file type cndb first $frame_num last $frame_num
set molid [molinfo top]
mol delrep 0 $molid

mol representation VDW 8.0 12.0
mol color Beta
mol selection "all"
mol material Opaque
mol addrep $molid

display projection orthographic
display depthcue off
color Display Background white
axes location off
display resetview

puts "Rendering view 3 (beads only, large)..."
render TachyonInternal ${output_base}_view3.tga

# Setup 3: Loops only
mol delete $molid
mol new $input_file type cndb first $frame_num last $frame_num
set molid [molinfo top]
mol delrep 0 $molid

mol representation Bonds 3.0 12.0
mol color ColorID 3
mol selection "all"
mol material Opaque
mol addrep $molid

display projection orthographic
display depthcue off
color Display Background white
axes location off
display resetview

puts "Rendering view 4 (loops only)..."
render TachyonInternal ${output_base}_view4.tga

puts ""
puts "Done! Created 4 views:"
puts "  ${output_base}_view1.tga - Beads + loops (beta coloring)"
puts "  ${output_base}_view2.tga - Beads + loops (occupancy coloring)"
puts "  ${output_base}_view3.tga - Beads only (large)"
puts "  ${output_base}_view4.tga - Loops only"
puts ""
puts "Convert to PNG with: convert image.tga image.png"

exit
