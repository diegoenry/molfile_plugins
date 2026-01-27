# Chromatin Tube Visualization for CNDB Files
# Self-contained script to draw smooth tubes through chromatin beads
# Based on VMD's NewCartoon B-spline approach
#
# Usage in VMD:
#   source chromatin_tube.tcl
#   draw_chromatin_tube top
#   draw_chromatin_tube top 10.0 16 8 blue Glossy

package require Tcl 8.5

###############################################################################
# Vector Math Utilities
###############################################################################

proc vadd {v1 v2} {
    lassign $v1 x1 y1 z1
    lassign $v2 x2 y2 z2
    return [list [expr {$x1 + $x2}] [expr {$y1 + $y2}] [expr {$z1 + $z2}]]
}

proc vsub {v1 v2} {
    lassign $v1 x1 y1 z1
    lassign $v2 x2 y2 z2
    return [list [expr {$x1 - $x2}] [expr {$y1 - $y2}] [expr {$z1 - $z2}]]
}

proc vscale {s v} {
    lassign $v x y z
    return [list [expr {$s * $x}] [expr {$s * $y}] [expr {$s * $z}]]
}

proc vdot {v1 v2} {
    lassign $v1 x1 y1 z1
    lassign $v2 x2 y2 z2
    return [expr {$x1*$x2 + $y1*$y2 + $z1*$z2}]
}

proc vcross {v1 v2} {
    lassign $v1 x1 y1 z1
    lassign $v2 x2 y2 z2
    return [list \
        [expr {$y1*$z2 - $z1*$y2}] \
        [expr {$z1*$x2 - $x1*$z2}] \
        [expr {$x1*$y2 - $y1*$x2}]]
}

proc vlength {v} {
    lassign $v x y z
    return [expr {sqrt($x*$x + $y*$y + $z*$z)}]
}

proc vnorm {v} {
    set len [vlength $v]
    if {$len < 1e-10} {
        return [list 0.0 0.0 0.0]
    }
    return [vscale [expr {1.0 / $len}] $v]
}

###############################################################################
# B-Spline Functions
###############################################################################

proc create_bspline_basis {} {
    return [list \
        [list [expr {-1.0/6.0}] [expr {3.0/6.0}] [expr {-3.0/6.0}] [expr {1.0/6.0}]] \
        [list [expr {3.0/6.0}] [expr {-6.0/6.0}] [expr {3.0/6.0}] [expr {0.0/6.0}]] \
        [list [expr {-3.0/6.0}] [expr {0.0/6.0}] [expr {3.0/6.0}] [expr {0.0/6.0}]] \
        [list [expr {1.0/6.0}] [expr {4.0/6.0}] [expr {1.0/6.0}] [expr {0.0/6.0}]] \
    ]
}

proc make_spline_Q_matrix {basis p0 p1 p2 p3} {
    set q [list]

    lassign $p0 x0 y0 z0
    lassign $p1 x1 y1 z1
    lassign $p2 x2 y2 z2
    lassign $p3 x3 y3 z3

    for {set i 0} {$i < 4} {incr i} {
        set basis_row [lindex $basis $i]

        set a [expr {[lindex $basis_row 0]*$x0 + [lindex $basis_row 1]*$x1 + \
                     [lindex $basis_row 2]*$x2 + [lindex $basis_row 3]*$x3}]
        set b [expr {[lindex $basis_row 0]*$y0 + [lindex $basis_row 1]*$y1 + \
                     [lindex $basis_row 2]*$y2 + [lindex $basis_row 3]*$y3}]
        set c [expr {[lindex $basis_row 0]*$z0 + [lindex $basis_row 1]*$z1 + \
                     [lindex $basis_row 2]*$z2 + [lindex $basis_row 3]*$z3}]

        lappend q [list $a $b $c]
    }

    return $q
}

proc make_spline_interpolation {q t} {
    lassign [lindex $q 0] q0x q0y q0z
    lassign [lindex $q 1] q1x q1y q1z
    lassign [lindex $q 2] q2x q2y q2z
    lassign [lindex $q 3] q3x q3y q3z

    set x [expr {$t * ($t * ($t * $q0x + $q1x) + $q2x) + $q3x}]
    set y [expr {$t * ($t * ($t * $q0y + $q1y) + $q2y) + $q3y}]
    set z [expr {$t * ($t * ($t * $q0z + $q1z) + $q2z) + $q3z}]

    return [list $x $y $z]
}

proc make_spline_tangent {q t} {
    lassign [lindex $q 0] q0x q0y q0z
    lassign [lindex $q 1] q1x q1y q1z
    lassign [lindex $q 2] q2x q2y q2z

    set t2 [expr {$t * $t}]

    set x [expr {3.0*$q0x*$t2 + 2.0*$q1x*$t + $q2x}]
    set y [expr {3.0*$q0y*$t2 + 2.0*$q1y*$t + $q2y}]
    set z [expr {3.0*$q0z*$t2 + 2.0*$q1z*$t + $q2z}]

    return [vnorm [list $x $y $z]]
}

proc generate_spline_points {atom_coords segments_per_span} {
    set n [llength $atom_coords]
    if {$n < 2} {
        error "Need at least 2 points to generate spline"
    }

    set spline_points [list]
    set tangents [list]
    set basis [create_bspline_basis]

    # Extend coordinates (duplicate endpoints)
    set extended_coords [list]
    lappend extended_coords [lindex $atom_coords 0]
    foreach coord $atom_coords {
        lappend extended_coords $coord
    }
    lappend extended_coords [lindex $atom_coords end]

    set splinedivs [expr {$segments_per_span + 2}]
    set inv_splinedivs [expr {1.0 / double($splinedivs)}]

    # Generate points for each span
    set n_extended [llength $extended_coords]
    for {set i 0} {$i < $n_extended - 3} {incr i} {
        set p0 [lindex $extended_coords $i]
        set p1 [lindex $extended_coords [expr {$i + 1}]]
        set p2 [lindex $extended_coords [expr {$i + 2}]]
        set p3 [lindex $extended_coords [expr {$i + 3}]]

        set q [make_spline_Q_matrix $basis $p0 $p1 $p2 $p3]

        for {set j 0} {$j < $splinedivs} {incr j} {
            set t [expr {$j * $inv_splinedivs}]
            set pt [make_spline_interpolation $q $t]
            set tan [make_spline_tangent $q $t]

            lappend spline_points $pt
            lappend tangents $tan
        }
    }

    # Add final point
    lappend spline_points [lindex $atom_coords end]
    set p0 [lindex $extended_coords [expr {$n_extended - 4}]]
    set p1 [lindex $extended_coords [expr {$n_extended - 3}]]
    set p2 [lindex $extended_coords [expr {$n_extended - 2}]]
    set p3 [lindex $extended_coords [expr {$n_extended - 1}]]
    set q [make_spline_Q_matrix $basis $p0 $p1 $p2 $p3]
    lappend tangents [make_spline_tangent $q 1.0]

    return [list $spline_points $tangents]
}

###############################################################################
# Rotation Minimizing Frames
###############################################################################

proc compute_rmf_frames {points tangents} {
    set n [llength $points]
    if {$n < 2} {
        error "Need at least 2 points for frame calculation"
    }

    set normals [list]
    set binormals [list]

    # Initialize first frame
    set t0 [lindex $tangents 0]

    if {abs([lindex $t0 0]) < 0.9} {
        set ref [list 1.0 0.0 0.0]
    } else {
        set ref [list 0.0 1.0 0.0]
    }

    set n0 [vnorm [vcross $t0 $ref]]
    set b0 [vnorm [vcross $t0 $n0]]

    lappend normals $n0
    lappend binormals $b0

    # Propagate frames
    for {set i 1} {$i < $n} {incr i} {
        set p_prev [lindex $points [expr {$i - 1}]]
        set p_curr [lindex $points $i]
        set t_prev [lindex $tangents [expr {$i - 1}]]
        set t_curr [lindex $tangents $i]
        set n_prev [lindex $normals [expr {$i - 1}]]
        set b_prev [lindex $binormals [expr {$i - 1}]]

        set v1 [vsub $p_curr $p_prev]
        set c1 [vdot $v1 $v1]

        if {$c1 < 1e-10} {
            lappend normals $n_prev
            lappend binormals $b_prev
            continue
        }

        set rL [vsub $n_prev [vscale [expr {2.0 / $c1 * [vdot $v1 $n_prev]}] $v1]]
        set tL [vsub $t_prev [vscale [expr {2.0 / $c1 * [vdot $v1 $t_prev]}] $v1]]

        set v2 [vsub $t_curr $tL]
        set c2 [vdot $v2 $v2]

        if {$c2 < 1e-10} {
            set n_curr $rL
        } else {
            set n_curr [vsub $rL [vscale [expr {2.0 / $c2 * [vdot $v2 $rL]}] $v2]]
        }

        set n_curr [vnorm $n_curr]
        set b_curr [vnorm [vcross $t_curr $n_curr]]

        lappend normals $n_curr
        lappend binormals $b_curr
    }

    return [list $normals $binormals]
}

###############################################################################
# Tube Drawing
###############################################################################

proc draw_tube_geometry {mol_id points normals binormals radius radial_segments color material {benchmark_var ""}} {
    set n_points [llength $points]
    set n_segments $radial_segments

    set PI 3.14159265359
    set frac_angle [expr {2.0 * $PI / double($n_segments)}]

    # Generate circles
    if {$benchmark_var ne ""} {
        upvar $benchmark_var bench
        set t0 [clock microseconds]
    }

    set circles [list]
    set circle_normals [list]

    for {set i 0} {$i < $n_points} {incr i} {
        set center [lindex $points $i]
        set normal [lindex $normals $i]
        set binormal [lindex $binormals $i]

        set circle [list]
        set cnorms [list]

        for {set j 0} {$j < $n_segments} {incr j} {
            set angle [expr {$j * $frac_angle}]
            set cos_a [expr {cos($angle)}]
            set sin_a [expr {sin($angle)}]

            set offset [vadd [vscale [expr {$radius * $cos_a}] $normal] \
                            [vscale [expr {$radius * $sin_a}] $binormal]]
            set vertex [vadd $center $offset]
            lappend circle $vertex

            set vnorm [vadd [vscale $cos_a $normal] [vscale $sin_a $binormal]]
            set vnorm [vnorm $vnorm]
            lappend cnorms $vnorm
        }

        lappend circles $circle
        lappend circle_normals $cnorms
    }

    if {$benchmark_var ne ""} {
        set t1 [clock microseconds]
        set bench(geometry_generation) [expr {($t1 - $t0) / 1000.0}]
    }

    # Draw triangles
    if {$benchmark_var ne ""} {
        set t0 [clock microseconds]
    }

    graphics $mol_id color $color
    graphics $mol_id material $material

    set n_triangles 0
    for {set i 0} {$i < $n_points - 1} {incr i} {
        set circle1 [lindex $circles $i]
        set circle2 [lindex $circles [expr {$i + 1}]]
        set norms1 [lindex $circle_normals $i]
        set norms2 [lindex $circle_normals [expr {$i + 1}]]

        for {set j 0} {$j < $n_segments} {incr j} {
            set j_next [expr {($j + 1) % $n_segments}]

            set v1 [lindex $circle1 $j]
            set v2 [lindex $circle1 $j_next]
            set v3 [lindex $circle2 $j_next]
            set v4 [lindex $circle2 $j]

            set n1 [lindex $norms1 $j]
            set n2 [lindex $norms1 $j_next]
            set n3 [lindex $norms2 $j_next]
            set n4 [lindex $norms2 $j]

            graphics $mol_id trinorm $v1 $v2 $v3 $n1 $n2 $n3
            graphics $mol_id trinorm $v1 $v3 $v4 $n1 $n3 $n4
            incr n_triangles 2
        }
    }

    if {$benchmark_var ne ""} {
        set t1 [clock microseconds]
        set bench(graphics_drawing) [expr {($t1 - $t0) / 1000.0}]
    }

    return $n_triangles
}

# Draw tube with colors varying by bead type
proc draw_tube_geometry_colored {mol_id points normals binormals radius radial_segments material natoms types resolution {benchmark_var ""}} {
    set n_points [llength $points]
    set n_segments $radial_segments

    # Define color mapping
    array set type_colors {
        A1 red
        A2 orange
        B1 blue
        B3 cyan
        NA gray
    }

    set PI 3.14159265359
    set frac_angle [expr {2.0 * $PI / double($n_segments)}]

    # Generate circles (same as before)
    if {$benchmark_var ne ""} {
        upvar $benchmark_var bench
        set t0 [clock microseconds]
    }

    set circles [list]
    set circle_normals [list]

    for {set i 0} {$i < $n_points} {incr i} {
        set center [lindex $points $i]
        set normal [lindex $normals $i]
        set binormal [lindex $binormals $i]

        set circle [list]
        set cnorms [list]

        for {set j 0} {$j < $n_segments} {incr j} {
            set angle [expr {$j * $frac_angle}]
            set cos_a [expr {cos($angle)}]
            set sin_a [expr {sin($angle)}]

            set offset [vadd [vscale [expr {$radius * $cos_a}] $normal] \
                            [vscale [expr {$radius * $sin_a}] $binormal]]
            set vertex [vadd $center $offset]
            lappend circle $vertex

            set vnorm [vadd [vscale $cos_a $normal] [vscale $sin_a $binormal]]
            set vnorm [vnorm $vnorm]
            lappend cnorms $vnorm
        }

        lappend circles $circle
        lappend circle_normals $cnorms
    }

    if {$benchmark_var ne ""} {
        set t1 [clock microseconds]
        set bench(geometry_generation) [expr {($t1 - $t0) / 1000.0}]
    }

    # Calculate spline points per bead span
    # splinedivs = resolution + 2 (as per generate_spline_points)
    set splinedivs [expr {$resolution + 2}]

    # Map each spline segment to its originating bead
    # The spline has (natoms-1) spans, each generating splinedivs points
    # Plus 1 final point

    graphics $mol_id material $material

    # Draw triangles with colors based on which bead span they belong to
    if {$benchmark_var ne ""} {
        set t0 [clock microseconds]
    }

    set n_triangles 0

    for {set i 0} {$i < $n_points - 1} {incr i} {
        set circle1 [lindex $circles $i]
        set circle2 [lindex $circles [expr {$i + 1}]]
        set norms1 [lindex $circle_normals $i]
        set norms2 [lindex $circle_normals [expr {$i + 1}]]

        # Determine which bead this segment belongs to
        # Each bead span creates splinedivs spline points
        # The first point duplicates endpoint, so we have (natoms+1) extended points
        # which create (natoms-1+2) = natoms+1 spans, but only natoms-1 original spans
        set bead_idx [expr {$i / $splinedivs}]

        # Clamp to valid range
        if {$bead_idx >= $natoms - 1} {
            set bead_idx [expr {$natoms - 2}]
        }
        if {$bead_idx < 0} {
            set bead_idx 0
        }

        # Get color for this bead
        set bead_type [lindex $types $bead_idx]
        if {[info exists type_colors($bead_type)]} {
            set seg_color $type_colors($bead_type)
        } else {
            set seg_color white
        }

        graphics $mol_id color $seg_color

        for {set j 0} {$j < $n_segments} {incr j} {
            set j_next [expr {($j + 1) % $n_segments}]

            set v1 [lindex $circle1 $j]
            set v2 [lindex $circle1 $j_next]
            set v3 [lindex $circle2 $j_next]
            set v4 [lindex $circle2 $j]

            set n1 [lindex $norms1 $j]
            set n2 [lindex $norms1 $j_next]
            set n3 [lindex $norms2 $j_next]
            set n4 [lindex $norms2 $j]

            graphics $mol_id trinorm $v1 $v2 $v3 $n1 $n2 $n3
            graphics $mol_id trinorm $v1 $v3 $v4 $n1 $n3 $n4
            incr n_triangles 2
        }
    }

    if {$benchmark_var ne ""} {
        set t1 [clock microseconds]
        set bench(graphics_drawing) [expr {($t1 - $t0) / 1000.0}]
    }

    return $n_triangles
}

###############################################################################
# Main Chromatin Tube Functions
###############################################################################

proc draw_chromatin_tube {mol_id {radius 10.0} {radial_segments 12} {resolution 6} {color "blue"} {material "Opaque"} {color_by_type 0}} {

    puts "\n=========================================="
    puts "Chromatin Tube Visualization"
    puts "=========================================="

    # Overall timing
    set t_total_start [clock microseconds]

    # Data preparation
    set t0 [clock microseconds]
    set sel [atomselect $mol_id "all"]
    set natoms [$sel num]
    set coords [$sel get {x y z}]

    # Get types if coloring by type
    if {$color_by_type} {
        set types [$sel get type]
        puts "Color mode: By compartment type"
    } else {
        set types [list]
        puts "Color mode: Single color ($color)"
    }
    $sel delete
    set t1 [clock microseconds]
    set time_data_prep [expr {($t1 - $t0) / 1000.0}]

    puts "Beads: $natoms"
    puts "Radius: $radius A"
    puts "Resolution: $resolution"

    if {$natoms < 2} {
        puts "ERROR: Need at least 2 beads"
        return
    }

    # B-spline generation
    puts "\nGenerating B-spline..."
    set t0 [clock microseconds]
    lassign [generate_spline_points $coords $resolution] spline_points tangents
    set t1 [clock microseconds]
    set time_spline [expr {($t1 - $t0) / 1000.0}]
    puts "Spline points: [llength $spline_points]"

    # RMF computation
    puts "Computing rotation minimizing frames..."
    set t0 [clock microseconds]
    lassign [compute_rmf_frames $spline_points $tangents] normals binormals
    set t1 [clock microseconds]
    set time_rmf [expr {($t1 - $t0) / 1000.0}]

    # Drawing
    puts "Drawing tube..."

    array set bench_times {}

    if {$color_by_type} {
        # Draw with varying colors based on bead types
        set n_tri [draw_tube_geometry_colored $mol_id $spline_points $normals $binormals \
                                               $radius $radial_segments $material \
                                               $natoms $types $resolution bench_times]
    } else {
        # Draw with single color
        set n_tri [draw_tube_geometry $mol_id $spline_points $normals $binormals \
                                       $radius $radial_segments $color $material bench_times]
    }

    set t_total_end [clock microseconds]
    set time_total [expr {($t_total_end - $t_total_start) / 1000.0}]

    # Print results
    puts "\n=========================================="
    puts "BENCHMARK RESULTS"
    puts "=========================================="
    puts [format "Triangles generated: %d" $n_tri]
    puts ""
    puts "Timing breakdown (milliseconds):"
    puts [format "  Data preparation:     %8.2f ms" $time_data_prep]
    puts [format "  B-spline generation:  %8.2f ms" $time_spline]
    puts [format "  RMF computation:      %8.2f ms" $time_rmf]
    if {[info exists bench_times(geometry_generation)]} {
        puts [format "  Geometry generation:  %8.2f ms" $bench_times(geometry_generation)]
    }
    if {[info exists bench_times(graphics_drawing)]} {
        puts [format "  Graphics drawing:     %8.2f ms" $bench_times(graphics_drawing)]
    }
    puts [format "  --------------------------------"]
    puts [format "  TOTAL TIME:           %8.2f ms" $time_total]
    puts ""

    # Calculate percentages
    puts "Time distribution:"
    set pct_data [expr {100.0 * $time_data_prep / $time_total}]
    set pct_spline [expr {100.0 * $time_spline / $time_total}]
    set pct_rmf [expr {100.0 * $time_rmf / $time_total}]
    puts [format "  Data preparation:     %5.1f%%" $pct_data]
    puts [format "  B-spline generation:  %5.1f%%" $pct_spline]
    puts [format "  RMF computation:      %5.1f%%" $pct_rmf]
    if {[info exists bench_times(geometry_generation)]} {
        set pct_geom [expr {100.0 * $bench_times(geometry_generation) / $time_total}]
        puts [format "  Geometry generation:  %5.1f%%" $pct_geom]
    }
    if {[info exists bench_times(graphics_drawing)]} {
        set pct_draw [expr {100.0 * $bench_times(graphics_drawing) / $time_total}]
        puts [format "  Graphics drawing:     %5.1f%%" $pct_draw]
    }

    puts "=========================================="
    puts "Complete!"
    puts "=========================================="
}

proc draw_chromatin_tube_by_type {mol_id {radius 10.0} {radial_segments 10} {resolution 4} {material "Opaque"}} {

    puts "\n=========================================="
    puts "Chromatin Tubes by Compartment Type"
    puts "=========================================="

    array set type_colors {
        A1 red
        A2 orange
        B1 blue
        B3 cyan
        NA gray
    }

    # Process each compartment type separately
    foreach type {A1 A2 B1 B3 NA} {
        set sel [atomselect $mol_id "type $type"]
        set natoms [$sel num]

        if {$natoms < 2} {
            $sel delete
            continue
        }

        set coords [$sel get {x y z}]
        set resids [$sel get resid]
        $sel delete

        puts "\nType $type: $natoms beads"

        # Group consecutive beads
        set segments [list]
        set current_segment [list [lindex $coords 0]]
        set prev_resid [lindex $resids 0]

        for {set i 1} {$i < $natoms} {incr i} {
            set curr_resid [lindex $resids $i]

            # If consecutive residues, add to segment
            if {$curr_resid == $prev_resid + 1} {
                lappend current_segment [lindex $coords $i]
            } else {
                # Save segment if it has enough points
                if {[llength $current_segment] >= 2} {
                    lappend segments $current_segment
                }
                set current_segment [list [lindex $coords $i]]
            }
            set prev_resid $curr_resid
        }

        # Add last segment
        if {[llength $current_segment] >= 2} {
            lappend segments $current_segment
        }

        # Draw each segment
        set color $type_colors($type)
        foreach segment $segments {
            if {[llength $segment] >= 2} {
                lassign [generate_spline_points $segment $resolution] pts tans
                lassign [compute_rmf_frames $pts $tans] norms binorms
                draw_tube_geometry $mol_id $pts $norms $binorms \
                                   $radius $radial_segments $color $material
            }
        }
    }

    puts "=========================================="
    puts "All compartments drawn!"
    puts "=========================================="
}

puts "=========================================="
puts "Chromatin Tube Drawer"
puts "=========================================="
puts "Commands:"
puts "  draw_chromatin_tube mol_id ?radius? ?segments? ?resolution? ?color? ?material? ?color_by_type?"
puts "  draw_chromatin_tube_by_type mol_id ?radius? ?segments? ?resolution? ?material?"
puts ""
puts "Default Parameters:"
puts "  radius: 10.0 A"
puts "  radial_segments: 12 (draw_chromatin_tube), 10 (by_type)"
puts "  resolution: 6 (draw_chromatin_tube), 4 (by_type)"
puts "  color_by_type: 0 = single color, 1 = color by compartment type"
puts ""
puts "Examples:"
puts "  # Single color tube (defaults: r=10, seg=12, res=6)"
puts "  draw_chromatin_tube top"
puts ""
puts "  # Custom parameters"
puts "  draw_chromatin_tube top 12.0 16 8 red Glossy"
puts ""
puts "  # Tube colored by compartment type (rainbow effect)"
puts "  draw_chromatin_tube top 10.0 12 6 blue Opaque 1"
puts ""
puts "  # Separate tubes for each compartment type"
puts "  draw_chromatin_tube_by_type top"
puts ""
puts "Clear: graphics top delete all"
puts "=========================================="
