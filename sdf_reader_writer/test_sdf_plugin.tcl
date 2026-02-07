# Test script for SDF/MOL V2000 plugin
# Usage: vmd -dispdev text -e test_sdf_plugin.tcl

vmdinfo version
puts "Loading SDF plugin..."

# Get plugin directory
set plugindir [file dirname [info script]]
set pluginpath "$plugindir/sdfplugin.so"

if {[file exists $pluginpath]} {
    puts "Found plugin at: $pluginpath"
    vmd_plugin_scandirectory $plugindir *.so
} else {
    puts "ERROR: Plugin not found at $pluginpath"
    puts "Please compile the plugin first with 'make'"
    quit
}

# Check if plugin was loaded
puts "\nChecking for SDF plugin..."
set found 0
for {set i 0} {$i < [molfile_plugin count]} {incr i} {
    set plugin_name [molfile_plugin name $i]
    if {[string match "*sdf*" $plugin_name]} {
        puts "  -> Found: $plugin_name"
        set found 1
    }
}

if {!$found} {
    puts "ERROR: SDF plugin not found after loading"
    quit
}

# Test loading the example SDF file
set testfile "$plugindir/example.sdf"
if {[file exists $testfile]} {
    puts "\n========================================="
    puts "Testing with file: $testfile"
    puts "=========================================\n"

    mol new $testfile type sdf waitfor all

    set molid [molinfo top]
    set natoms [molinfo $molid get numatoms]
    set nframes [molinfo $molid get numframes]

    puts "Molecule loaded successfully!"
    puts "  Number of atoms: $natoms"
    puts "  Number of frames: $nframes"

    # Verify atom count (ethanol = 9 atoms)
    if {$natoms != 9} {
        puts "FAIL: Expected 9 atoms, got $natoms"
    } else {
        puts "PASS: Atom count correct (9)"
    }

    # Check element names
    set sel [atomselect $molid "index 0"]
    set name [$sel get name]
    puts "\nFirst atom name: $name"
    if {$name eq "C"} {
        puts "PASS: First atom is carbon"
    } else {
        puts "FAIL: Expected C, got $name"
    }
    $sel delete

    set sel [atomselect $molid "index 2"]
    set name [$sel get name]
    puts "Third atom name: $name"
    if {$name eq "O"} {
        puts "PASS: Third atom is oxygen"
    } else {
        puts "FAIL: Expected O, got $name"
    }
    $sel delete

    # Check bonds
    set nbonds [molinfo $molid get numbonds]
    puts "\nNumber of bonds: $nbonds"
    if {$nbonds == 8} {
        puts "PASS: Bond count correct (8)"
    } else {
        puts "FAIL: Expected 8 bonds, got $nbonds"
    }

    # Check coordinates of first atom
    set sel [atomselect $molid "index 0"]
    set x [$sel get x]
    set y [$sel get y]
    set z [$sel get z]
    puts "\nFirst atom coords: x=$x y=$y z=$z"
    if {abs($x - (-0.7516)) < 0.001} {
        puts "PASS: X coordinate correct"
    } else {
        puts "FAIL: Expected x=-0.7516, got $x"
    }
    $sel delete

    puts "\n========================================="
    puts "Plugin test completed!"
    puts "=========================================\n"

} else {
    puts "ERROR: Test file not found: $testfile"
}

# Uncomment to quit automatically
# quit
