# Chromatin Tube Drawer for CNDB files
# Draws a smooth tube through chromatin beads in sequence
# Based on VMD's NewCartoon representation using B-splines
#
# Usage:
#   source draw_chromatin_tube.tcl
#   draw_chromatin_tube top ?radius? ?radial_segments? ?resolution? ?color? ?material?

package require Tcl 8.5

# Source the worm tube library
set scriptdir [file dirname [info script]]
if {[file exists "$scriptdir/../Worms/worm_tube.tcl"]} {
    source "$scriptdir/../Worms/worm_tube.tcl"
} else {
    puts "ERROR: Cannot find worm_tube.tcl"
    puts "Please ensure worm_tube.tcl is in ../Worms/ directory"
    return
}

# Draw chromatin tube through all beads in sequence
# Treats all chromatin beads as a connected polymer chain
proc draw_chromatin_tube {mol_id {radius 10.0} {radial_segments 16} {resolution 10} {color "blue"} {material "Opaque"}} {

    puts "\n=========================================="
    puts "Chromatin Tube Visualization"
    puts "=========================================="
    puts "Molecule: $mol_id"
    puts "Radius: $radius A"
    puts "Radial segments: $radial_segments"
    puts "Resolution: $resolution"
    puts "Color: $color"
    puts "Material: $material"
    puts ""

    # Get all chromatin beads in sequence (by resid order)
    set sel [atomselect $mol_id "all"]
    set natoms [$sel num]
    puts "Total beads: $natoms"

    if {$natoms < 2} {
        puts "ERROR: Need at least 2 beads to draw tube"
        $sel delete
        return
    }

    # Get coordinates sorted by resid (bead sequence)
    set coords [$sel get {x y z}]
    $sel delete

    puts "Drawing tube through $natoms beads..."
    puts ""

    # Use the worm tube function with B-spline
    # Note: We pass coordinates directly since beads are already in sequence order
    set atom_coords $coords
    set use_bspline 1  # Use B-spline for smoothness

    # Generate spline points
    puts "Generating B-spline..."
    lassign [generate_spline_points $atom_coords $resolution $use_bspline] spline_points tangents
    set n_spline [llength $spline_points]
    puts "Generated $n_spline spline points"

    # Compute rotation minimizing frames
    puts "Computing rotation minimizing frames..."
    lassign [compute_rmf_frames $spline_points $tangents] normals binormals

    # Draw tube geometry
    puts "Drawing tube geometry..."
    set n_triangles [draw_tube_geometry $mol_id $spline_points $tangents $normals $binormals \
                                        $radius $radial_segments $color $material]

    puts "Drew $n_triangles triangles"
    puts "=========================================="
    puts "Chromatin tube complete!"
    puts "=========================================="
}

# Draw tubes colored by compartment type
proc draw_chromatin_tubes_by_type {mol_id {radius 10.0} {radial_segments 16} {resolution 10} {material "Opaque"}} {

    puts "\n=========================================="
    puts "Chromatin Tubes by Compartment Type"
    puts "=========================================="

    # Define colors for each type
    array set type_colors {
        A1 red
        A2 orange
        B1 blue
        B3 cyan
        NA gray
    }

    foreach type {A1 A2 B1 B3 NA} {
        set sel [atomselect $mol_id "type $type"]
        set natoms [$sel num]

        if {$natoms < 2} {
            puts "Skipping $type: only $natoms beads"
            $sel delete
            continue
        }

        puts "\nDrawing tube for $type ($natoms beads)..."

        # Get coordinates for this type
        set coords [$sel get {x y z}]
        $sel delete

        # For each segment of consecutive beads of the same type,
        # we need to check if they're actually consecutive in sequence
        # For now, we'll draw all beads of this type as segments

        # Use worm tube with this type's color
        set color $type_colors($type)

        # Generate spline
        lassign [generate_spline_points $coords $resolution 1] spline_points tangents

        # Compute frames
        lassign [compute_rmf_frames $spline_points $tangents] normals binormals

        # Draw tube
        set n_triangles [draw_tube_geometry $mol_id $spline_points $tangents $normals $binormals \
                                            $radius $radial_segments $color $material]

        puts "  $type: $n_triangles triangles"
    }

    puts "=========================================="
    puts "All compartment tubes complete!"
    puts "=========================================="
}

# Draw tube with variable radius based on compartment type
proc draw_chromatin_tube_variable_radius {mol_id {base_radius 10.0} {radial_segments 16} {resolution 10} {material "Opaque"}} {

    puts "\n=========================================="
    puts "Chromatin Tube with Variable Radius"
    puts "=========================================="

    # Get all beads
    set sel [atomselect $mol_id "all"]
    set natoms [$sel num]
    set types [$sel get type]
    set coords [$sel get {x y z}]
    $sel delete

    puts "Total beads: $natoms"

    # Draw segments with different radii
    # A compartments: thicker (more open/active)
    # B compartments: thinner (more compact/inactive)
    array set radius_scale {
        A1 1.2
        A2 1.1
        B1 0.9
        B3 0.8
        NA 0.7
    }

    array set type_colors {
        A1 red
        A2 orange
        B1 blue
        B3 cyan
        NA gray
    }

    # Group consecutive beads and draw segments
    set start_idx 0
    set current_type [lindex $types 0]

    for {set i 1} {$i <= $natoms} {incr i} {
        set next_type ""
        if {$i < $natoms} {
            set next_type [lindex $types $i]
        }

        # If type changed or end of sequence, draw segment
        if {$next_type ne $current_type || $i == $natoms} {
            set end_idx $i
            set segment_length [expr {$end_idx - $start_idx}]

            if {$segment_length >= 2} {
                # Extract coordinates for this segment
                set segment_coords [lrange $coords $start_idx [expr {$end_idx - 1}]]

                # Calculate radius for this type
                set seg_radius [expr {$base_radius * $radius_scale($current_type)}]
                set color $type_colors($current_type)

                puts "Drawing $current_type segment: beads $start_idx-[expr {$end_idx-1}] (radius: [format %.1f $seg_radius] A)"

                # Generate and draw
                lassign [generate_spline_points $segment_coords $resolution 1] spline_points tangents
                lassign [compute_rmf_frames $spline_points $tangents] normals binormals
                draw_tube_geometry $mol_id $spline_points $tangents $normals $binormals \
                                   $seg_radius $radial_segments $color $material
            }

            set start_idx $i
            set current_type $next_type
        }
    }

    puts "=========================================="
    puts "Variable radius tube complete!"
    puts "=========================================="
}

puts "=========================================="
puts "Chromatin Tube Drawer Loaded"
puts "=========================================="
puts "Functions:"
puts "  draw_chromatin_tube mol_id ?radius? ?radial_segments? ?resolution? ?color? ?material?"
puts "  draw_chromatin_tubes_by_type mol_id ?radius? ?radial_segments? ?resolution? ?material?"
puts "  draw_chromatin_tube_variable_radius mol_id ?base_radius? ?radial_segments? ?resolution? ?material?"
puts ""
puts "Examples:"
puts "  # Basic tube through all beads"
puts "  draw_chromatin_tube top 10.0 16 10 blue Opaque"
puts ""
puts "  # Tubes colored by compartment type"
puts "  draw_chromatin_tubes_by_type top 10.0 16 10 Opaque"
puts ""
puts "  # Variable radius based on compartment type"
puts "  draw_chromatin_tube_variable_radius top 10.0 16 10 Opaque"
puts ""
puts "Clear graphics: graphics top delete all"
puts "=========================================="
