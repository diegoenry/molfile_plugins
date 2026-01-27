#!/usr/bin/env vmd
# Example analysis script for CNDB files
# Demonstrates how to perform common chromatin dynamics analyses

# Load plugin
set scriptdir [file dirname [info script]]
if {[file exists "$scriptdir/cndbplugin.so"]} {
    vmd_plugin_scandirectory $scriptdir *.so
    puts "Loaded CNDB plugin"
}

# Configuration
if {$argc >= 1} {
    set cndb_file [lindex $argv 0]
} else {
    set cndb_file "/Users/deb0054/github/NDB-Converters/Bcell_short.cndb"
}

if {$argc >= 2} {
    set nframes [lindex $argv 1]
} else {
    set nframes 100
}

if {$argc >= 3} {
    set output_file [lindex $argv 2]
} else {
    set output_file "analysis_results.dat"
}

puts "=========================================="
puts "CNDB Chromatin Dynamics Analysis"
puts "=========================================="
puts "File: $cndb_file"
puts "Frames: $nframes"
puts ""

# Load structure and trajectory
puts "Loading data..."
mol new $cndb_file type cndb first 0 last [expr $nframes - 1] waitfor all
set molid [molinfo top]

# Get basic information
set natoms [molinfo $molid get numatoms]
set nframes_loaded [molinfo $molid get numframes]

puts "Loaded $natoms atoms, $nframes_loaded frames"
puts ""

# Analysis 1: Radius of Gyration vs Time
puts "Computing radius of gyration..."
proc compute_rg {molid frame} {
    set sel [atomselect $molid "all" frame $frame]
    set coords [$sel get {x y z}]
    set n [llength $coords]

    # Center of mass
    set cx 0.0; set cy 0.0; set cz 0.0
    foreach coord $coords {
        set cx [expr $cx + [lindex $coord 0]]
        set cy [expr $cy + [lindex $coord 1]]
        set cz [expr $cz + [lindex $coord 2]]
    }
    set cx [expr $cx / $n]
    set cy [expr $cy / $n]
    set cz [expr $cz / $n]

    # Rg
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

set rg_list {}
for {set i 0} {$i < $nframes_loaded} {incr i} {
    set rg [compute_rg $molid $i]
    lappend rg_list $rg
    if {$i % 20 == 0} {
        puts "  Frame $i: Rg = [format %.3f $rg] A"
    }
}

# Statistics
set rg_sum 0.0
foreach rg $rg_list {
    set rg_sum [expr $rg_sum + $rg]
}
set rg_mean [expr $rg_sum / [llength $rg_list]]
puts "  Mean Rg: [format %.3f $rg_mean] A"
puts ""

# Analysis 2: End-to-End Distance
puts "Computing end-to-end distance..."
set ete_list {}
for {set i 0} {$i < $nframes_loaded} {incr i} {
    set sel1 [atomselect $molid "index 0" frame $i]
    set sel2 [atomselect $molid "index [expr $natoms - 1]" frame $i]
    set pos1 [lindex [$sel1 get {x y z}] 0]
    set pos2 [lindex [$sel2 get {x y z}] 0]
    set dist [veclength [vecsub $pos2 $pos1]]
    lappend ete_list $dist
    $sel1 delete
    $sel2 delete
    if {$i % 20 == 0} {
        puts "  Frame $i: d(0-[expr $natoms-1]) = [format %.3f $dist] A"
    }
}

set ete_sum 0.0
foreach ete $ete_list {
    set ete_sum [expr $ete_sum + $ete]
}
set ete_mean [expr $ete_sum / [llength $ete_list]]
puts "  Mean end-to-end distance: [format %.3f $ete_mean] A"
puts ""

# Analysis 3: Contact Map from Loop Constraints
puts "Analyzing chromatin loops..."
set sel [atomselect $molid "all"]
set bondlist [$sel getbonds]
$sel delete

# Count contacts
set contact_map {}
for {set i 0} {$i < $natoms} {incr i} {
    set bonds [lindex $bondlist $i]
    foreach j $bonds {
        if {$i < $j} {  # Avoid double counting
            lappend contact_map [list $i $j]
        }
    }
}

puts "  Total loop contacts: [llength $contact_map]"

# Find longest-range contact
set max_range 0
set max_pair {0 0}
foreach contact $contact_map {
    set i [lindex $contact 0]
    set j [lindex $contact 1]
    set range [expr abs($j - $i)]
    if {$range > $max_range} {
        set max_range $range
        set max_pair [list $i $j]
    }
}
puts "  Longest-range contact: $max_pair (range: $max_range beads)"
puts ""

# Analysis 4: Average Loop Distance
puts "Computing average loop distances..."
set loop_distances {}
for {set frame 0} {$frame < $nframes_loaded} {incr frame} {
    foreach contact $contact_map {
        set i [lindex $contact 0]
        set j [lindex $contact 1]
        set sel1 [atomselect $molid "index $i" frame $frame]
        set sel2 [atomselect $molid "index $j" frame $frame]
        set pos1 [lindex [$sel1 get {x y z}] 0]
        set pos2 [lindex [$sel2 get {x y z}] 0]
        set dist [veclength [vecsub $pos2 $pos1]]
        lappend loop_distances $dist
        $sel1 delete
        $sel2 delete
    }
}

set loop_sum 0.0
foreach dist $loop_distances {
    set loop_sum [expr $loop_sum + $dist]
}
set loop_mean [expr $loop_sum / [llength $loop_distances]]
puts "  Mean loop distance: [format %.3f $loop_mean] A"
puts "  (averaged over all loops and all frames)"
puts ""

# Analysis 5: Genomic Coverage
puts "Analyzing genomic coverage..."
set sel [atomselect $molid "all"]
set beta_vals [$sel get beta]
set occ_vals [$sel get occupancy]
$sel delete

set min_pos [lindex [lsort -real $beta_vals] 0]
set max_pos [lindex [lsort -real $occ_vals] end]
set coverage [expr $max_pos - $min_pos]
puts "  Genomic range: [format %.0f $min_pos] - [format %.0f $max_pos] bp"
puts "  Coverage: [format %.1f [expr $coverage / 1e6]] Mb"
puts "  Resolution: [format %.1f [expr $coverage / $natoms / 1e3]] kb/bead"
puts ""

# Save results to file
puts "Saving results to $output_file..."
set fout [open $output_file w]
puts $fout "# CNDB Chromatin Dynamics Analysis Results"
puts $fout "# File: $cndb_file"
puts $fout "# Date: [clock format [clock seconds]]"
puts $fout "#"
puts $fout "# Atoms: $natoms"
puts $fout "# Frames analyzed: $nframes_loaded"
puts $fout "# Genomic range: [format %.0f $min_pos] - [format %.0f $max_pos] bp"
puts $fout "#"
puts $fout "# Frame\tRg(A)\tEte(A)"
for {set i 0} {$i < $nframes_loaded} {incr i} {
    puts $fout "$i\t[format %.4f [lindex $rg_list $i]]\t[format %.4f [lindex $ete_list $i]]"
}
close $fout

puts "Results saved!"
puts ""
puts "=========================================="
puts "Summary Statistics"
puts "=========================================="
puts "  Atoms: $natoms beads"
puts "  Genomic coverage: [format %.1f [expr $coverage / 1e6]] Mb"
puts "  Mean Rg: [format %.3f $rg_mean] A"
puts "  Mean end-to-end: [format %.3f $ete_mean] A"
puts "  Loop contacts: [llength $contact_map]"
puts "  Mean loop distance: [format %.3f $loop_mean] A"
puts "=========================================="
puts ""

# Quit
quit
