# Benchmark script for chromatin tube visualization
# Tests different resolutions and parameters to identify bottlenecks

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

puts "\n############################################################"
puts "CHROMATIN TUBE BENCHMARKING SUITE"
puts "############################################################"
puts ""
puts "This script tests different parameter combinations to"
puts "identify computational bottlenecks in tube visualization."
puts ""

# Test 1: Default parameters (baseline)
puts "\n============================================================"
puts "TEST 1: Default parameters (single color)"
puts "Parameters: radius=10, radial_segments=12, resolution=6"
puts "============================================================"
graphics $molid delete all
draw_chromatin_tube $molid

# Test 2: Colored by type
puts "\n\n============================================================"
puts "TEST 2: Rainbow tube (colored by type)"
puts "Parameters: radius=10, radial_segments=12, resolution=6, color_by_type=1"
puts "============================================================"
graphics $molid delete all
draw_chromatin_tube $molid 10.0 12 6 blue Opaque 1

# Test 3: High resolution
puts "\n\n============================================================"
puts "TEST 3: High resolution (single color)"
puts "Parameters: radius=10, radial_segments=16, resolution=8"
puts "============================================================"
graphics $molid delete all
draw_chromatin_tube $molid 10.0 16 8 blue Opaque 0

# Test 4: Low resolution
puts "\n\n============================================================"
puts "TEST 4: Low resolution (single color)"
puts "Parameters: radius=10, radial_segments=8, resolution=4"
puts "============================================================"
graphics $molid delete all
draw_chromatin_tube $molid 10.0 8 4 blue Opaque 0

# Test 5: Very high spline resolution (stress test)
puts "\n\n============================================================"
puts "TEST 5: Very high spline resolution (stress test)"
puts "Parameters: radius=10, radial_segments=12, resolution=12"
puts "============================================================"
graphics $molid delete all
draw_chromatin_tube $molid 10.0 12 12 blue Opaque 0

# Test 6: Very high radial segments (stress test)
puts "\n\n============================================================"
puts "TEST 6: Very high radial segments (stress test)"
puts "Parameters: radius=10, radial_segments=24, resolution=6"
puts "============================================================"
graphics $molid delete all
draw_chromatin_tube $molid 10.0 24 6 blue Opaque 0

puts "\n\n############################################################"
puts "BENCHMARK COMPLETE"
puts "############################################################"
puts ""
puts "ANALYSIS:"
puts "Compare the timing percentages across different tests to"
puts "identify which operations scale poorly with parameters."
puts ""
puts "- Spline generation scales with 'resolution' parameter"
puts "- Geometry generation scales with both parameters"
puts "- Graphics drawing scales with total triangle count"
puts ""
puts "To clear all graphics: graphics top delete all"
puts "############################################################"
