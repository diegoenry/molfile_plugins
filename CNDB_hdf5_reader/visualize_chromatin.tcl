# Visualization script for OpenMiChroM CNDB files
# Usage: vmd -e visualize_chromatin.tcl -args your_file.cndb [nframes]
#
# Example: vmd -e visualize_chromatin.tcl -args Bcell_short.cndb 100

# Load CNDB plugin from current directory
set scriptdir [file dirname [info script]]
if {[file exists "$scriptdir/cndbplugin.so"]} {
    vmd_plugin_scandirectory $scriptdir *.so
    puts "Loaded CNDB plugin from $scriptdir"
} else {
    puts "Warning: cndbplugin.so not found in $scriptdir"
    puts "Assuming plugin is already installed in VMD"
}

# Get command line arguments
if {$argc >= 1} {
    set cndb_file [lindex $argv 0]
} else {
    puts "Usage: vmd -e visualize_chromatin.tcl -args <cndb_file> \[nframes\]"
    puts "Example: vmd -e visualize_chromatin.tcl -args Bcell_short.cndb 100"
    exit
}

# Number of frames to load (default: 100)
if {$argc >= 2} {
    set nframes [lindex $argv 1]
} else {
    set nframes 100
}

puts "========================================"
puts "OpenMiChroM Chromatin Visualization"
puts "========================================"
puts "File: $cndb_file"
puts "Loading first $nframes frames..."
puts ""

# Load the CNDB file
mol new $cndb_file type cndb first 0 last [expr $nframes - 1] waitfor all

# Get molecule info
set molid [molinfo top]
set natoms [molinfo $molid get numatoms]
set nframes_loaded [molinfo $molid get numframes]

puts "Loaded molecule:"
puts "  Atoms (beads): $natoms"
puts "  Frames: $nframes_loaded"
puts ""

# Delete default representation
mol delrep 0 $molid

# Representation 1: Chromatin beads colored by genomic position
puts "Creating representations..."
mol representation VDW 5.0 12.0
mol color Beta
mol selection "all"
mol material Opaque
mol addrep $molid
puts "  1. Chromatin beads (colored by genomic position)"

# Representation 2: Chromatin loops (bonds)
mol representation Bonds 2.0 12.0
mol color ColorID 3
mol selection "all"
mol material Opaque
mol addrep $molid
puts "  2. Chromatin loops (cyan bonds)"

# Representation 3: Backbone connections (optional - connects sequential beads)
# Uncomment to show polymer backbone
# mol representation Bonds 1.0 12.0
# mol color ColorID 1
# mol selection "all"
# mol material Opaque
# mol addrep $molid
# puts "  3. Polymer backbone (red bonds)"

# Set display properties
display projection orthographic
display depthcue off
color Display Background white
axes location off

# Center molecule
display resetview

puts ""
puts "Visualization controls:"
puts "  - Mouse: Rotate (left), Translate (middle), Scale (right)"
puts "  - Animation: Use VMD animation controls or type 'animate goto N'"
puts "  - To save image: render TachyonInternal image.tga"
puts ""

# Print genomic range
set sel [atomselect $molid "all"]
set beta_vals [$sel get beta]
set occ_vals [$sel get occupancy]
set min_pos [lindex [lsort -real $beta_vals] 0]
set max_pos [lindex [lsort -real $occ_vals] end]
puts "Genomic coordinate range: [format %.0f $min_pos] - [format %.0f $max_pos] bp"
$sel delete

# Compute radius of gyration for first frame
proc compute_rg {molid frameid} {
    set sel [atomselect $molid "all" frame $frameid]
    set coords [$sel get {x y z}]
    set n [llength $coords]

    # Compute center of mass
    set cx 0.0
    set cy 0.0
    set cz 0.0
    foreach coord $coords {
        set cx [expr $cx + [lindex $coord 0]]
        set cy [expr $cy + [lindex $coord 1]]
        set cz [expr $cz + [lindex $coord 2]]
    }
    set cx [expr $cx / $n]
    set cy [expr $cy / $n]
    set cz [expr $cz / $n]

    # Compute Rg
    set sum_sq 0.0
    foreach coord $coords {
        set dx [expr [lindex $coord 0] - $cx]
        set dy [expr [lindex $coord 1] - $cy]
        set dz [expr [lindex $coord 2] - $cz]
        set sum_sq [expr $sum_sq + ($dx*$dx + $dy*$dy + $dz*$dz)]
    }
    $sel delete
    return [expr sqrt($sum_sq / $n)]
}

puts ""
puts "Computing radius of gyration..."
set rg0 [compute_rg $molid 0]
puts "  Frame 0: [format %.2f $rg0] Angstroms"

if {$nframes_loaded > 1} {
    set rg_last [compute_rg $molid [expr $nframes_loaded - 1]]
    puts "  Frame [expr $nframes_loaded - 1]: [format %.2f $rg_last] Angstroms"
}

puts ""
puts "========================================"
puts "Visualization ready!"
puts "========================================"
puts ""

# Optional: Start animation
# animate goto 0
# animate forward
