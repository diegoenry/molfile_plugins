# Test script for CNDB plugin
# Usage: vmd -dispdev text -e test_cndb_plugin.tcl

# Load the plugin
vmdinfo version
puts "Loading CNDB plugin..."

# Get plugin directory
set plugindir [file dirname [info script]]
set pluginpath "$plugindir/cndbplugin.so"

if {[file exists $pluginpath]} {
    puts "Found plugin at: $pluginpath"
    # Load the plugin library
    vmd_plugin_scandirectory $plugindir *.so
} else {
    puts "ERROR: Plugin not found at $pluginpath"
    puts "Please compile the plugin first with 'make'"
    quit
}

# Check if plugin was loaded
puts "\nChecking for CNDB plugin..."
set found 0
foreach type {molfile} {
    for {set i 0} {$i < [eval "${type}_plugin count"]} {incr i} {
        set plugin_name [eval "${type}_plugin name $i"]
        if {[string match "*cndb*" $plugin_name]} {
            puts "  -> Found: $plugin_name"
            set found 1
        }
    }
}

if {!$found} {
    puts "ERROR: CNDB plugin not found after loading"
    puts "Available molfile plugins:"
    for {set i 0} {$i < [molfile_plugin count]} {incr i} {
        puts "  - [molfile_plugin name $i]"
    }
    quit
}

# Test loading the example file
set testfile "/Users/deb0054/github/NDB-Converters/Bcell_short.cndb"
if {[file exists $testfile]} {
    puts "\n========================================="
    puts "Testing with file: $testfile"
    puts "=========================================\n"

    # Load only the structure first (no trajectory)
    mol new $testfile type cndb first 0 last 0 waitfor all

    # Get molecule info
    set molid [molinfo top]
    set natoms [molinfo $molid get numatoms]
    set nframes [molinfo $molid get numframes]

    puts "\nMolecule loaded successfully!"
    puts "  Number of atoms: $natoms"
    puts "  Number of frames: $nframes"

    # Get first atom info
    set sel [atomselect $molid "index 0"]
    puts "\nFirst atom properties:"
    puts "  Name: [$sel get name]"
    puts "  Type: [$sel get type]"
    puts "  Resname: [$sel get resname]"
    puts "  Resid: [$sel get resid]"
    puts "  Chain: [$sel get chain]"
    puts "  Segid: [$sel get segid]"
    puts "  Beta (genomic start): [$sel get beta]"
    puts "  Occupancy (genomic end): [$sel get occupancy]"
    $sel delete

    # Get coordinates from first and last frame
    animate goto 0
    set sel [atomselect $molid "index 0"]
    set coords0 [$sel get {x y z}]
    puts "\nFirst atom coordinates (frame 0): $coords0"
    $sel delete

    if {$nframes > 1} {
        animate goto [expr $nframes - 1]
        set sel [atomselect $molid "index 0"]
        set coordsN [$sel get {x y z}]
        puts "First atom coordinates (frame [expr $nframes - 1]): $coordsN"
        $sel delete
    }

    # Check bonds (loop constraints)
    set nbonds [molinfo $molid get numbonds]
    puts "\nNumber of bonds (loop constraints): $nbonds"
    if {$nbonds > 0} {
        puts "First 5 bonds:"
        set sel [atomselect $molid "all"]
        set bondlist [$sel getbonds]
        for {set i 0} {$i < 5 && $i < [llength $bondlist]} {incr i} {
            set bonds [lindex $bondlist $i]
            if {[llength $bonds] > 0} {
                puts "  Atom $i bonded to: $bonds"
            }
        }
        $sel delete
    }

    puts "\n========================================="
    puts "Plugin test completed successfully!"
    puts "=========================================\n"

} else {
    puts "ERROR: Test file not found: $testfile"
}

# Uncomment to quit automatically
# quit
