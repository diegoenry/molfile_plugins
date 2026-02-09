# test_mmcif_plugin.tcl
#
# VMD test script for the mmCIF reader plugin.
# Run inside VMD:  source test_mmcif_plugin.tcl
#
# Requirements:
#   - mmcifplugin.so in the current directory
#   - 1ubq.cif downloaded (the script will try to fetch it if missing)

proc test_section {name} {
    puts "\n=========================================="
    puts "TEST: $name"
    puts "=========================================="
}

proc test_assert {desc condition} {
    if {$condition} {
        puts "  PASS: $desc"
    } else {
        puts "  FAIL: $desc"
    }
}

# Load the plugin
test_section "Plugin Registration"

set plugindir [pwd]
vmd_plugin_scandirectory $plugindir

# Check if plugin is registered
set found 0
foreach p [molfile list] {
    if {[lindex $p 0] eq "mmcif"} {
        set found 1
        puts "  Found plugin: $p"
    }
}
test_assert "mmcif plugin registered" $found

# Download test file if not present
if {![file exists 1ubq.cif]} {
    puts "\nDownloading 1UBQ.cif from RCSB PDB..."
    set url "https://files.rcsb.org/download/1UBQ.cif"
    if {[catch {
        package require http
        package require tls
        ::http::register https 443 ::tls::socket
        set tok [::http::geturl $url]
        set fp [open 1ubq.cif w]
        puts -nonewline $fp [::http::data $tok]
        close $fp
        ::http::cleanup $tok
        puts "  Downloaded 1ubq.cif"
    } err]} {
        puts "  Could not download: $err"
        puts "  Please download manually:"
        puts "    curl -o 1ubq.cif https://files.rcsb.org/download/1UBQ.cif"
        puts "  Then re-run this script."
        return
    }
}

# Test 1: Basic loading
test_section "Basic Loading (1UBQ)"

mol new 1ubq.cif type mmcif waitfor all

set molid [molinfo top]
set natoms [molinfo $molid get numatoms]
puts "  Loaded molecule $molid with $natoms atoms"

# 1UBQ has 660 atoms (no hydrogens, single model)
test_assert "atom count is 660" [expr {$natoms == 660}]

# Check chain ID
set sel [atomselect $molid "all"]
set chains [lsort -unique [$sel get chain]]
$sel delete
test_assert "single chain A" [expr {$chains eq "A"}]

# Check residue names
set sel [atomselect $molid "resname MET and name CA"]
set count [$sel num]
$sel delete
test_assert "has MET residues with CA" [expr {$count >= 1}]

# Check element assignments
set sel [atomselect $molid "element C"]
set ncarbon [$sel num]
$sel delete
puts "  Carbon atoms: $ncarbon"
test_assert "has carbon atoms (element assignment works)" [expr {$ncarbon > 0}]

set sel [atomselect $molid "element N"]
set nnitrogen [$sel num]
$sel delete
puts "  Nitrogen atoms: $nnitrogen"
test_assert "has nitrogen atoms" [expr {$nnitrogen > 0}]

# Check residue IDs
set sel [atomselect $molid "resid 1"]
set n_res1 [$sel num]
$sel delete
test_assert "residue 1 exists" [expr {$n_res1 > 0}]

set sel [atomselect $molid "resid 76"]
set n_res76 [$sel num]
$sel delete
test_assert "residue 76 exists (last residue)" [expr {$n_res76 > 0}]

# Check coordinates are reasonable (not all zeros)
set sel [atomselect $molid "index 0"]
set coords [lindex [$sel get {x y z}] 0]
$sel delete
set x [lindex $coords 0]
set y [lindex $coords 1]
set z [lindex $coords 2]
puts "  First atom coords: $x $y $z"
test_assert "coordinates are non-zero" [expr {$x != 0.0 || $y != 0.0 || $z != 0.0}]

mol delete $molid

# Test 2: Chain filtering
test_section "Chain Filtering"

# 1UBQ only has chain A, so filtering for A should give same count
mol new "1ubq.cif?chain=A" type mmcif waitfor all
set molid [molinfo top]
set natoms_filtered [molinfo $molid get numatoms]
puts "  Filtered (chain=A): $natoms_filtered atoms"
test_assert "chain A filter gives 660 atoms" [expr {$natoms_filtered == 660}]
mol delete $molid

# Filtering for non-existent chain should give 0 atoms (or fail to load)
if {[catch {
    mol new "1ubq.cif?chain=Z" type mmcif waitfor all
    set molid [molinfo top]
    set natoms_z [molinfo $molid get numatoms]
    puts "  Filtered (chain=Z): $natoms_z atoms"
    mol delete $molid
} err]} {
    puts "  Chain Z filter correctly failed to load: $err"
}

# Test 3: Unit cell
test_section "Unit Cell"

mol new 1ubq.cif type mmcif waitfor all
set molid [molinfo top]
set a [molinfo $molid get a]
set b [molinfo $molid get b]
set c [molinfo $molid get c]
puts "  Unit cell: a=$a b=$b c=$c"
test_assert "unit cell A > 0" [expr {$a > 0}]
test_assert "unit cell B > 0" [expr {$b > 0}]
test_assert "unit cell C > 0" [expr {$c > 0}]
mol delete $molid

puts "\n=========================================="
puts "All tests complete."
puts "=========================================="
