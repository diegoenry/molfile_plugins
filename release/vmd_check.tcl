# Load release plugins into a real VMD and read a file with each one.
#
#   VMDDIR=/path/to/vmd /path/to/vmd/vmd_ARCH -dispdev text \
#       -e release/vmd_check.tcl -args PLUGIN_FOLDER WORK_FOLDER
#
# PLUGIN_FOLDER holds the .so files (for example $W/out/MACOSXARM64/molfile or
# an extracted archive's molfile/); WORK_FOLDER is make_release.sh's work
# folder ($W), for its test files. Prints RESULT lines and quits.
set dir [lindex $argv 0]
set W [lindex $argv 1]
set T $W/testfiles
set S $W/src
set failed 0

foreach so {pdbxplugin sdfplugin cndbplugin h5mdplugin_gromacs meshhdf5plugin} {
  if {[catch {plugin dlopen [file join $dir $so.so]} r] || $r < 1} {
    puts "RESULT FAILED dlopen $so: $r"
    incr failed
  }
}
foreach t {pdbx mmcif sdf cndb h5md_gromacs h5mesh} {
  catch {unset info}
  if {[catch {plugin info "mol file reader" $t info}] || ![info exists info(majorversion)]} {
    puts "RESULT FAILED plugin $t is not registered"
    incr failed
  } else {
    puts "RESULT plugin $t $info(majorversion).$info(minorversion)"
  }
}

# Each file is opened by its extension, unless a type is given.
proc check {label file args} {
  global failed
  if {[catch {eval mol new [list $file] $args waitfor all} m]} {
    puts "RESULT FAILED $label: $m"
    incr failed
    return
  }
  set n [molinfo $m get numatoms]
  set line "natoms=$n frames=[molinfo $m get numframes]"
  if {$n > 0} {
    set sel [atomselect $m all]
    set nb 0
    foreach b [$sel getbonds] { incr nb [llength $b] }
    append line " bonds=[expr {$nb / 2}] center=[format {%.3f %.3f %.3f} {*}[measure center $sel]]"
    $sel delete
  }
  append line " graphics=[llength [graphics $m list]]"
  puts "RESULT $label $line"
  mol delete $m
}

check "pdbx .cif" $S/pdbx/tests/data/chains.cif
check "pdbx .bcif" $S/pdbx/tests/bcif/chains.bcif
check "pdbx .cif.gz" $T/chains.cif.gz type pdbx
check "pdbx ?chain=A" "$S/pdbx/tests/data/chains.cif?chain=A" type pdbx
check "sdf" $S/sdf_reader_writer/example.sdf
check "h5md_gromacs" $S/gromacs_h5md/equil.h5md
set f [lindex [lsort [glob -nocomplain $T/*.cndb]] 0]
if {$f ne ""} { check "cndb (5 frames)" $f first 0 last 4 } else { puts "RESULT skipped cndb: no test file" }
set f [lindex [lsort [glob -nocomplain $T/*.h5mesh]] 0]
if {$f ne ""} { check "h5mesh" $f } else { puts "RESULT skipped h5mesh: no test file" }
puts "RESULT $failed failed"
quit
